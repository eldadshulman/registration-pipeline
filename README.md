# HEST-VALIS: H&E to Xenium single-cell registration + QC

Align an H&E whole-slide image onto its matched Xenium spatial dataset, at single-cell
accuracy, and certify the alignment with quantitative QC. Each slide automatically gets the
registration setting (micro non-rigid refinement, or not) that aligns it best.

The output is a **warped H&E image in the Xenium coordinate frame** plus a **per-slide QC
report**, so every transcript / cell can be placed on the right piece of tissue.

---

## What it does

1. **Register** the H&E onto the Xenium **DAPI** (the fixed reference) with VALIS using the
   HEST `register_dapi_he` recipe (H&E color-deconvolved, rigid + non-rigid), so the H&E ends
   up in the Xenium/molecular frame.
2. **Optionally refine** with `register_micro` (local non-rigid). Both variants are produced
   from a single registration.
3. **QC** each variant: warp the H&E nuclei (StarDist) into the Xenium frame and compare to the
   Xenium nuclei. Four checks (below).
4. **Select** micro vs no-micro **per slide** by a simple rule, so the cohort is a mix.
5. **Warp the full H&E image** with the chosen setting into the Xenium frame -> registered
   OME-TIFF.

Direction matters: DAPI is fixed and the **H&E moves onto it**, so all outputs share the
Xenium coordinate system. (Some pipelines do the reverse and only transfer annotations onto
cells; this one produces a real registered H&E image.)

## The four QC checks

| check | plain meaning | metric |
|-------|---------------|--------|
| nucleus coincidence | do the same nuclei line up? | median offset (um), target < 10 (one cell) |
| density correlation | do dense regions match dense regions? | Pearson `density_r` on 50 um bins |
| tissue occupancy | do cells sit on tissue, not glass? | fraction of cells on the tissue mask |
| negative control | does the QC break when we break the alignment? | density-r collapse under a +/-100 um shift |

## The per-slide selection rule (`hest_valis/select.py`)

- Primary metric is **nucleus-coincidence median um** -- lower wins.
- **Tie** (medians within 0.15 um): take the higher `density_r`.
- **Over-fit guard**: if the lower-um winner's `density_r` is > 0.10 below the loser's AND its
  um advantage is < 0.5 um, take the other protocol instead.

A negative `density_r` is a red flag (alignment wrong at the field level) even if the median
um passes -- quarantine and re-register such slides.

---

## Quickstart

```bash
cd hest_valis_pipeline
cp examples/config.json config.json        # edit pixel_um, output_dir
cp examples/samples.csv  samples.csv        # one row per slide (paths below)
# edit the env paths at the top of slurm/*.sbatch  (see env/setup.md)

# 1) QC + per-slide selection (array over samples; 0-indexed, header skipped)
sbatch --array=0-$(( $(tail -n +2 samples.csv | wc -l) - 1 )) slurm/qc_array.sbatch

# 2) aggregate the decisions
python run_select.py --samples samples.csv --config config.json
#    -> output/per_slide_decision.csv , output/wsi_manifest.csv

# 3) warp the chosen-protocol H&E image per slide (slow, CPU; see env/setup.md)
sbatch --array=0-$(( $(tail -n +2 output/wsi_manifest.csv | wc -l) - 1 )) slurm/wsi_array.sbatch
```

`samples.csv` columns:

| column | what |
|--------|------|
| `sample_id` | unique name; becomes the per-sample output folder |
| `he_path` | H&E whole-slide image (`.svs` / `.ome.tiff`), the moving image |
| `dapi_path` | Xenium DAPI `morphology_focus/ch0000_dapi.ome.tif`, the fixed reference |
| `xenium_cells` | Xenium `cells.parquet` (centroids in microns) for QC |

## Outputs (under `output/<sample_id>/`)

```
he_nuclei.npy              StarDist H&E nuclei (H&E pixels)
he_nuclei_nomicro.npy      H&E nuclei warped into the Xenium frame (no micro)
he_nuclei_micro.npy        ... with micro refinement (absent if micro failed)
qc.json                    both variants' metrics + the chosen protocol + the rule fired
registered/aligned_fullres_HE.ome.tiff   the warped H&E in the Xenium frame
```
Cohort level (under `output/`): `per_slide_decision.csv`, `wsi_manifest.csv`.

## Run without SLURM (one slide, interactively)

```bash
export PYTHONPATH=$PWD
stardist_env/bin/python  run_segment.py  --samples samples.csv --config config.json --sample SLIDE_A
valis_env/bin/python     run_register.py --samples samples.csv --config config.json --sample SLIDE_A
qc_env/bin/python        run_qc.py       --samples samples.csv --config config.json --sample SLIDE_A
valis_env/bin/python     run_wsi.py      --samples samples.csv --config config.json --sample SLIDE_A
```

## Use the library directly

```python
from hest_valis import registration, segment, concordance, select, xenium
reg = registration.register_slide(he, dapi, work_dir, micro=True)
warped = registration.warp_points(reg, segment.segment_he(he))   # nuclei -> Xenium frame
m = concordance.compute_qc(warped, xenium.load_xenium_nuclei(cells, 0.2125), 0.2125)
```

## Layout

```
hest_valis/        registration, segment, concordance, select, xenium, config
run_segment.py     step 1  (StarDist env)   H&E nuclei
run_register.py    step 2  (valis env)      register + warp nuclei, both protocols
run_qc.py          step 3  (QC env)         metrics + per-slide selection
run_wsi.py                 (valis env)      warp the chosen H&E image -> OME-TIFF
run_select.py      aggregate decisions -> decision table + WSI manifest
slurm/             SLURM array wrappers
env/setup.md       build the envs + the REQUIRED serial-read patch
examples/          config.json + samples.csv templates
```

## Gotchas

- **Image warp deadlock** -> apply the serial-read patch in `env/setup.md`. Without it,
  `run_wsi.py` hangs at "COLLECTING RESULTS" (thread-unsafe BioFormats reader, not a memory
  problem). The image warp is then slow (~15 tiles/s); give `wsi_array.sbatch` a generous time.
- **`aligned_fullres_HE`** -> valis hard-codes this moving-slide name; the pipeline symlinks
  your H&E to it automatically. Don't rename.
- **Registrar pickle is not reloadable** -> register and warp in the same process (the scripts
  do this). To get both QC variants cheaply, `run_register.py` warps nuclei before and after
  `register_micro` from one registration.
- **Pixel size** -> set `pixel_um` to your DAPI um/pixel (Xenium is 0.2125).

## Credits

This pipeline stands on two pieces of work and would not exist without them:

- **VALIS** -- the whole-slide image registration engine used here.
  Gatenbee et al., "Virtual Alignment of pathoLogy Image Series for multi-gene analysis,"
  *Nature Communications* (2023). https://github.com/MathOnco/valis
- **HEST / Mahmood Lab** -- the `register_dapi_he` recipe and the HEST-1k spatial
  transcriptomics + histology resource that this registration approach is based on.
  Jaume et al., "HEST-1k: A Dataset for Spatial Transcriptomics and Histology Image Analysis,"
  *NeurIPS* (2024). Mahmood Lab, https://github.com/mahmoodlab/HEST

Nuclei are detected with **StarDist** (Schmidt et al., MICCAI 2018). Please cite VALIS, HEST,
and StarDist if you use this pipeline.
