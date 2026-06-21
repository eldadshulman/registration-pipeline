"""HEST-VALIS registration: align an H&E whole-slide image onto the Xenium DAPI frame.

This is the HEST `register_dapi_he` recipe with ONE deliberate deviation: the Xenium DAPI
is the fixed reference, so the H&E is warped onto the molecular (Xenium) coordinate frame.
That means everything downstream (warped H&E image, warped nuclei) lands in Xenium pixel
space, ready for transcript / cell assignment.

Recipe:
  - Valis(reference_img_f=DAPI, align_to_reference=True, check_for_reflections=False)
  - register(): rigid + non-rigid (OpticalFlow), H&E color-deconvolved (HEDeconvolution),
    both slides read with BioFormats (thread-unsafe single reader, see env/setup.md).
  - register_micro(): optional local non-rigid refinement at up to 10000 px.

Two outputs you can take from one registration (no need to register twice):
  - warp_points(): warp a set of (x, y) points (e.g. H&E nuclei) into the Xenium frame. Cheap.
  - warp_image(): warp the full H&E image into the Xenium frame and save an OME-TIFF. Expensive.

IMPORTANT (the aligned_fullres_HE name): valis serialises the registrar mid-run and the HEST
fork hard-codes the slide name 'aligned_fullres_HE'. So the H&E file passed to valis MUST be
named 'aligned_fullres_HE.<ext>'. register_slide() handles this by symlinking your H&E to that
name inside the per-sample work dir. Do not rename it away.
"""
import os

# valis_hest pulls in a JVM + BioFormats; it is imported lazily inside the functions that
# need it (register_slide / shutdown) so this module stays importable -- and unit-testable
# (e.g. ometiff_pages) -- without the heavy registration env.

MICRO_MAX_DIM_PX = 10000  # HEST default for register_micro


def _aligned_he_symlink(he_path, work_dir):
    """valis/HEST hard-code the moving-slide name 'aligned_fullres_HE'; give it that name."""
    ext = os.path.splitext(he_path)[1] or ".svs"
    link = os.path.join(work_dir, f"aligned_fullres_HE{ext}")
    if not (os.path.islink(link) or os.path.exists(link)):
        os.symlink(os.path.abspath(he_path), link)
    return link


def register_slide(he_path, dapi_path, work_dir, micro=False):
    """Register H&E onto DAPI. Returns an in-process valis registrar.

    he_path   : H&E whole-slide image (e.g. .svs / .ome.tiff), the MOVING image.
    dapi_path : Xenium DAPI (morphology_focus ch0), the FIXED reference.
    work_dir  : per-sample scratch dir (valis writes here).
    micro     : if True, also run register_micro() for local non-rigid refinement.

    The registrar pickle valis writes is NOT reliably reloadable, so keep this registrar
    object in-process and warp from it directly (warp_points / warp_image below).
    """
    from valis_hest import preprocessing, registration
    from valis_hest.slide_io import BioFormatsSlideReader
    os.makedirs(work_dir, exist_ok=True)
    reg_dir = os.path.join(work_dir, "valis_output")
    os.makedirs(reg_dir, exist_ok=True)
    he = _aligned_he_symlink(he_path, work_dir)

    reg = registration.Valis(
        "", reg_dir,
        img_list=[dapi_path, he],
        reference_img_f=dapi_path,         # DAPI is the fixed reference (deviation from HEST)
        align_to_reference=True,
        check_for_reflections=False,
    )
    reg.register(
        brightfield_processing_cls=preprocessing.HEDeconvolution,
        reader_dict={dapi_path: [BioFormatsSlideReader], he: [BioFormatsSlideReader]},
    )
    if micro:
        reg.register_micro(
            max_non_rigid_registration_dim_px=MICRO_MAX_DIM_PX,
            align_to_reference=True,
            brightfield_processing_cls=preprocessing.HEDeconvolution,
            reference_img_f=dapi_path,
        )
    return reg


def he_slide(reg):
    """Return the H&E (moving) slide object from a registrar."""
    return [s for nm, s in reg.slide_dict.items() if "dapi" not in nm.lower()][0]


def warp_points(reg, xy, non_rigid=True):
    """Warp (N, 2) points from H&E pixel space into the Xenium/DAPI frame.

    With micro registration applied, non_rigid=True uses the micro-refined field; without it,
    the register() OpticalFlow field. This avoids the expensive (and historically fragile)
    full-image warp, so it is the right tool for QC.
    """
    return he_slide(reg).warp_xy(xy, slide_level=0, pt_level=0, non_rigid=non_rigid, crop="reference")


def ometiff_pages(path):
    """Number of IFD pages in an OME-TIFF. 0 means a truncated/incomplete write: a job killed
    mid-write leaves a multi-GB file whose page table was never finalised (BigTIFF first-IFD
    offset still 0), so no reader can open it. Use this to tell a real WSI from a corpse."""
    import tifffile
    try:
        with tifffile.TiffFile(path) as tf:
            return len(tf.pages)
    except Exception:
        return 0


def warp_image(reg, out_dir, level=0, non_rigid=True, compression="deflate"):
    """Warp the full H&E image into the Xenium frame and save an OME-TIFF in out_dir.

    Requires the serial-read patch in env/setup.md or valis deadlocks at the COLLECTING step.

    Crash-safe write: valis streams a pyramidal OME-TIFF, so a job killed mid-write (preemption
    / walltime / OOM) leaves a multi-GB file with an empty page table -- unreadable, yet present,
    so a naive "skip if exists" would never regenerate it. To avoid that, warp into a sibling
    .tmp dir, verify every output actually opens (pages > 0), then atomically move into place.
    An interrupted run therefore leaves NO file at the canonical path, so the next run re-warps
    instead of skipping a corpse.
    """
    import glob
    os.makedirs(out_dir, exist_ok=True)
    tmp_dir = out_dir.rstrip("/") + ".tmp"
    if os.path.isdir(tmp_dir):
        for f in glob.glob(os.path.join(tmp_dir, "*")):
            try:
                os.remove(f)
            except OSError:
                pass
    os.makedirs(tmp_dir, exist_ok=True)

    reg.warp_and_save_slides(tmp_dir, level=level, non_rigid=non_rigid, crop="reference",
                             compression=compression)

    written = glob.glob(os.path.join(tmp_dir, "*.ome.tif*"))
    if not written:
        raise RuntimeError(f"warp_image: valis wrote no OME-TIFF into {tmp_dir}")
    for f in written:
        if ometiff_pages(f) == 0:
            raise RuntimeError(f"warp_image: {os.path.basename(f)} has 0 pages (truncated write); "
                               f"leaving {out_dir} empty so the next run re-warps")
    for f in written:                                   # all good -> publish atomically
        os.replace(f, os.path.join(out_dir, os.path.basename(f)))
    try:
        os.rmdir(tmp_dir)
    except OSError:
        pass
    return out_dir


def prerotate_he(he_path, rotation, out_path):
    """Losslessly rotate the H&E by a cardinal angle (0/90/180/270 clockwise) with pyvips and
    save a tiled pyramidal TIFF. Used to seed registration for a grossly mis-oriented slide
    (see coarse_align.cardinal_rotation). Returns (out_path, (width, height))."""
    import pyvips
    img = pyvips.Image.new_from_file(he_path, access="sequential")
    img = {0: img, 90: img.rot90(), 180: img.rot180(), 270: img.rot270()}[rotation]
    img.tiffsave(out_path, tile=True, pyramid=True, compression="jpeg", Q=92, bigtiff=True)
    return out_path, (img.width, img.height)


def shutdown():
    from valis_hest import registration
    try:
        registration.kill_jvm()
    except Exception:
        pass
