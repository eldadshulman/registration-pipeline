# Environment setup

This pipeline uses three environments. They can be three conda envs, or two (if your
StarDist and QC deps coexist). Only the **valis** env needs the patch below.

## 1. valis env (registration + image warp)

HEST ships a vendored fork of valis as the `valis_hest` package. Build an env with it plus a
working Java (BioFormats uses jpype/JVM):

```bash
conda create -n valis_hest python=3.10 openjdk maven -y
conda activate valis_hest
pip install valis-wsi            # or install HEST and use its bundled valis_hest
pip install jpype1 pyvips tifffile zarr pandas numpy scipy
```

Point `VALIS_PY` / `VALIS_ENV` in the SLURM wrappers at this env.

### REQUIRED PATCH: serial tile read (or the image warp deadlocks)

`valis_hest/slide_io.py` reads slide tiles with many threads sharing **one** BioFormats
reader. BioFormats readers are not thread-safe, so the threads deadlock at the "COLLECTING
RESULTS" step and `warp_and_save_slides` hangs forever (it is NOT a memory problem; more RAM
does not help). The fix is to serialise the tile read.

Find `get_tiles_parallel` in `slide_io.py` and change:

```python
        n_cpu = valtils.get_ncpus_available() - 1
```
to:
```python
        n_cpu = 1   # BioFormats reader is not thread-safe; serial read avoids the COLLECTING deadlock
```

Or apply the patch:

```bash
SLIDE_IO=$VALIS_ENV/lib/python3.10/site-packages/valis_hest/slide_io.py
cp "$SLIDE_IO" "$SLIDE_IO.orig"
sed -i 's/n_cpu = valtils.get_ncpus_available() - 1/n_cpu = 1  # serial read: BioFormats reader not thread-safe/' "$SLIDE_IO"
```

Tradeoff: the warp is slower (~15 tiles/s, so a 100k-tile slide reads in ~2 h), but reliable.
This only affects `run_wsi.py` (the full image warp). The QC path (`run_register.py`, which
warps nuclei *points*) never hits this and needs no patch. Give `wsi_array.sbatch` a generous
`--time` (default 14 h).

Note: the registrar also hard-codes the moving-slide name `aligned_fullres_HE`. The pipeline
handles this automatically by symlinking your H&E to that name in each sample's work dir
(`registration._aligned_he_symlink`) -- do not rename it away.

## 2. StarDist env (H&E nuclei)

```bash
conda create -n stardist python=3.10 -y
conda activate stardist
pip install stardist tensorflow csbdeep scikit-image tifffile zarr numpy
```

## 3. QC env (metrics + selection)

Any env with `numpy scipy pandas tifffile zarr`. The valis or stardist env usually already
satisfies this, so you may reuse one of them for `run_qc.py`.
