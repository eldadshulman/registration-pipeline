# HEST-VALIS: H&E to Xenium single-cell registration + QC

Align an H&E whole-slide image onto its matched Xenium spatial dataset, **targeting** single-cell
accuracy (a median nucleus offset below ~10 um, the size of one cell) and certifying it with
quantitative QC. Each slide automatically gets the registration setting (micro non-rigid
refinement, or not) that aligns it best.

The output is a **warped H&E image in the Xenium coordinate frame** plus a **per-slide QC
report**, so every transcript / cell can be placed on the right piece of tissue.

> Built on [VALIS](https://github.com/MathOnco/valis), the HEST
> [`register_dapi_he`](https://github.com/mahmoodlab/HEST) recipe, and
> [StarDist](https://github.com/stardist/stardist) (the `2D_versatile_he` model). **If you use
> this pipeline, please also cite the upstream tools** -- see [Credits](#credits) and
> `CITATION.cff`.

Why this is needed: 10x's own documentation notes the post-Xenium H&E is imaged on a *different*
microscope and is **not pre-registered** to the Xenium data
([Understanding Xenium Outputs](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs)),
so it must be registered before transcripts can be placed on H&E tissue.

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
5. **Warp the full H&E image** with the chosen setting into the Xenium frame -> a registered
   [OME-TIFF](https://ome-model.readthedocs.io/en/stable/ome-tiff/).
6. **Transfer per-cell annotations**: derive an H&E region map (tumor / stroma / background)
   and tag every Xenium cell with the region it falls in.

If a slide comes back with a **negative density-r** (a gross mis-orientation that VALIS could
not recover), a **coarse rotation/flip search runs automatically as a fallback** and rescues
it (see Self-healing below).

Direction matters: DAPI is fixed and the **H&E moves onto it**, so all outputs share the
Xenium coordinate system. This pipeline produces both a real registered H&E image *and*
per-cell annotations.

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

# 4) per-slide + cohort alignment QC report (CPU-only plotting; no re-registration)
sbatch --array=0-$(( $(tail -n +2 samples.csv | wc -l) - 1 )) slurm/report_array.sbatch
#    -> output/<sample>/report.pdf  and  output/cohort_report.pdf (the triage page)

# 5) acceptance gate + move/audit table (same qc.json + thresholds as the cohort report)
python run_provenance.py --samples samples.csv --config config.json
#    -> output/provenance.csv  (accepted vs manual-review, per-slide audit)
```

#### One-command cohort driver (`run_cohort.sh`)

To run all five steps end to end -- in order, **blocking on each array before the next** -- use
`run_cohort.sh`. It can't be a plain `&&` chain (each `sbatch` returns immediately, and the
`wsi_array` size isn't known until `run_select` writes `wsi_manifest.csv`), so it launches each
SLURM stage with `sbatch --wait` and sizes the wsi array right before submitting it.

```bash
module load slurm/slurm-compbio/23.11.10          # sbatch on PATH (Cedars)
SELECT_PY=qc_env/bin/python \
  ./run_cohort.sh --samples samples.csv --config config.json
#   --dry-run      print the plan + array sizes, submit nothing
#   --strict       abort the cohort if any array task fails (default: warn + continue)
#   --skip-wsi / --skip-report   drop those stages (e.g. QC-only cohort pass)
```

It blocks for as long as the warp takes (hours) -- run it under `tmux`/`nohup`, or submit it as a
tiny long-lived `defq` job (it just waits). Individual slide failures don't abort by default:
`run_select`/`run_provenance` aggregate whatever completed and report the misses. The array TASKS
read the `samples.csv`/`config.json` in `slurm/*.sbatch` (their EDIT-THESE block) -- keep the
`--samples`/`--config` here pointing at the same files. Covered by `tests/test_run_cohort_sh.py`
(stubbed `sbatch` + interpreter).

`samples.csv` columns:

| column | what |
|--------|------|
| `sample_id` | unique name; becomes the per-sample output folder |
| `he_path` | H&E whole-slide image (`.svs` / `.ome.tiff`), the moving image |
| `dapi_path` | Xenium DAPI image, the fixed reference. Raw 10x output is `morphology_focus/morphology_focus_0000.ome.tif` (channel 0 = DAPI); HEST renames it to `morphology_focus/ch0000_dapi.ome.tif`. Point at whichever you have. |
| `xenium_cells` | Xenium [`cells.parquet`](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs) (cell centroids in microns) for QC |
| `batch` *(optional)* | scanner-run key; `run_provenance.py` uses it for per-batch orientation-consistency checks. Absent -> all slides share batch `UNKNOWN`. |

## Outputs (under `output/<sample_id>/`)

```
he_nuclei.npy              StarDist H&E nuclei (H&E pixels)
he_nuclei_nomicro.npy      H&E nuclei warped into the Xenium frame (no micro)
he_nuclei_micro.npy        ... with micro refinement (absent if micro failed)
he_nuclei_coarse.npy       ... coarse-fallback alignment (only if a rescue was needed)
qc.json                    all variants' metrics + the chosen protocol + the rule fired
registered/aligned_fullres_HE.ome.tiff   the warped H&E in the Xenium frame
cell_labels.parquet        per-cell annotation: cell_id, x_um, y_um, he_region
region_overlay.png         tumor/stroma/background region map (QC)
report.pdf                 one-page alignment QC report (+ report.png for the review notebook)
```
Cohort level (under `output/`): `per_slide_decision.csv`, `wsi_manifest.csv`, `cohort_report.pdf`,
`provenance.csv`.

## One triage source: thresholds + provenance

The cohort report (`cohort_report.pdf`) and the provenance/move table (`provenance.csv`) read the
**same `qc.json`** and gate on the **same thresholds**, so they can never tell you different things
about the same slide.

- **One schema.** Every slide's outcome lives in `qc.json` (`decision.chosen` +
  `decision.sel_density_r` / `sel_median_um`, with the chosen metric under `metrics[<chosen>]`).
  A rescue is not a separate file: `run_rescue.py` writes `decision.chosen="rescued"` and
  `metrics["rescued"]` in the `compute_qc` shape, so the report and the provenance pick up rescued
  slides for free.
- **One set of thresholds** (in `config.json`, the single source -- `hest_valis/config.py`):

  | key | used by | meaning |
  |-----|---------|---------|
  | `density_r_accept` | report cohort line **and** provenance gate | accept if `density_r >= this` |
  | `median_um_accept` | report cohort line **and** provenance gate | accept if median offset `<= this` (um) |
  | `rescue_trigger_r` | `run_qc.py` | selected `density_r` below this -> attempt the coarse/orient rescue |
  | `rescue_delta_min` | provenance | a rescued slide whose `density_r` jump is below this is flagged for eyes |

  The report's cohort "good" line is literally `density_r_accept` / `median_um_accept`, the same
  numbers `provenance.gate()` accepts on -- change them in one place and both move together.

```bash
# acceptance gate + audit/move table (CPU; same qc.json + thresholds as the cohort report)
python run_provenance.py --samples samples.csv --config config.json
#    -> output/provenance.csv  + a printed triage summary (accepted vs manual-review)
```

`provenance.csv` is one row per slide: the gate result (`accepted`), the as-is-vs-rescued audit
(`pre_r/post_r`, `delta_r`, recovered orientation), `reason`, per-batch `orientation_outlier`
(needs an optional `batch` column in `samples.csv`), `small_delta_flag`, and the chosen artifacts
(`chosen_nuclei`, `registered_wsi`). A slide that fails the gate is flagged `manual`, never
auto-accepted.

## Optional flags (`run_qc.py`)

| flag | effect |
|------|--------|
| `--no-occupancy` | skip tissue-mask computation and the occupancy QC check (faster; use when DAPI quality is poor) |
| `--no-coarse-fallback` | disable the automatic coarse rotation/flip rescue for slides with negative density-r |

## Run without SLURM (one slide, interactively)

```bash
export PYTHONPATH=$PWD
stardist_env/bin/python  run_segment.py  --samples samples.csv --config config.json --sample SLIDE_A
valis_env/bin/python     run_register.py --samples samples.csv --config config.json --sample SLIDE_A
qc_env/bin/python        run_qc.py       --samples samples.csv --config config.json --sample SLIDE_A
valis_env/bin/python     run_wsi.py      --samples samples.csv --config config.json --sample SLIDE_A
```

### One-command single-slide wrapper (`run_slide.py`) -- demo / debugging

`run_slide.py` chains those steps for **one slide** -- register -> QC -> selection -> (optional)
WSI warp -- and gates the expensive warp on the pipeline's own decision. It does **not**
re-implement any QC threshold, selection, or rescue logic: it invokes the same entry points and
consumes the `qc.json` decision, importing `config.thresholds` / `provenance.gate` for the accept
terminology. **This is a local convenience for demos/debugging; the documented production path
remains the SLURM arrays above.**

```bash
export PYTHONPATH=$PWD
# per-stage envs differ (StarDist / valis / QC); each falls back to the current python if unset
export STARDIST_PY=stardist_env/bin/python  VALIS_PY=valis_env/bin/python  QC_PY=qc_env/bin/python

# register + QC + selection, then STOP and print what would happen (no warp):
python run_slide.py --samples samples.csv --config config.json --sample SLIDE_A
# add the expensive full-res warp -- only if selection marks the slide eligible:
python run_slide.py --samples samples.csv --config config.json --sample SLIDE_A --warp
python run_slide.py ... --dry-run          # print the plan + intended actions; run nothing
python run_slide.py ... --resume           # skip stages whose valid outputs already exist
python run_slide.py ... --force-stage qc   # re-run a stage even if its output exists (repeatable)
```

It prints a stage summary (status, selected protocol, QC values, decision reason, whether the
warp ran, output paths) and returns a documented exit code:

| code | status | meaning |
|---|---|---|
| 0 | `ELIGIBLE_FOR_WARP` | reached an eligible selection (and warped, if `--warp`) |
| 2 | `REVIEW_REQUIRED` | selected but below the accept gate, or an orientation rescue (warp via `run_rescue.py --warp-image`) |
| 3 | `QC_FAILED` | no usable registration decision |
| 1 | (stage failed) | a stage command exited non-zero -- execution stops |

The warp runs **only** with `--warp` **and** a `micro`/`nomicro` selection that passes the accept
gate; `coarse`/`rescued` slides are reported for `run_rescue.py --warp-image`. Orchestration is
covered by `tests/test_run_slide.py` (mocked stages).

#### Plain shell chain (`run_slide.sh`)

If you want the transparent version with **no decision logic** -- just `run_segment.py -> run_register.py
-> run_qc.py` chained with `&&`, each in its own env -- use `run_slide.sh`. It is exactly equivalent
to running the three commands by hand: it stops on the first failure and its exit status is the
failing stage's status (no remapped `2`/`3` codes, no `qc.json` read-back, no warp). Selection still
happens inside `run_qc.py`. Per-stage interpreters come from `$STARDIST_PY`/`$VALIS_PY`/`$QC_PY`.

```bash
STARDIST_PY=stardist_env/bin/python VALIS_PY=valis_env/bin/python QC_PY=qc_env/bin/python \
  ./run_slide.sh --samples samples.csv --config config.json --sample SLIDE_A
```

Covered by `tests/test_run_slide_sh.py` (stubbed interpreters).

## Use the library directly

```python
from hest_valis import registration, segment, concordance, select, xenium
reg = registration.register_slide(he, dapi, work_dir, micro=True)
warped = registration.warp_points(reg, segment.segment_he(he))   # nuclei -> Xenium frame
m = concordance.compute_qc(warped, xenium.load_xenium_nuclei(cells, 0.2125), 0.2125)
```

## Per-cell annotation transfer (`run_annotate.py`)

After registration, tag every Xenium cell with its H&E region. A region map (tumor / stroma /
background) is built from the registered-H&E nuclear density (2-component GMM: the higher-density
tissue cluster is tumor, the lower is stroma, empty bins are background), then each Xenium cell is
assigned the region it falls in.

```bash
python run_annotate.py --samples samples.csv --config config.json --sample SLIDE_A
# -> output/SLIDE_A/cell_labels.parquet   (cell_id, x_um, y_um, he_region)
#    output/SLIDE_A/region_overlay.png
```

This is the same idea as annotation-transfer pipelines that overlap each cell with an aligned
mask, except the mask is derived from H&E morphology. **When real pathologist masks are
available**, replace `annotate.region_map()` with a lookup into that mask; `assign_cells()` is
unchanged, and you get pathologist-grade per-cell labels.

## Alignment QC report (`run_report.py`)

A **non-interactive** per-slide diagnostic PDF (Agg backend, no prompts) so the whole cohort can
be eyeballed without opening each slide by hand. It only **reads** what the pipeline already
produced -- it never re-registers, and only reloads a whole-slide image when the warped raster
already exists. CPU-only (just plotting). The matching and binning are reused from
`concordance` (`mutual_nn_pairs`, `density_grids`), so the report and the QC numbers agree by
construction.

```bash
# per slide -> output/<sample>/report.pdf (+ report.png)
python run_report.py --samples samples.csv --config config.json --sample SLIDE_A
# cohort triage page -> output/cohort_report.pdf
python run_report.py --samples samples.csv --config config.json --cohort
# whole cohort on SLURM (CPU; last task also renders the cohort page):
sbatch --array=0-$(( $(tail -n +2 samples.csv | wc -l) - 1 )) slurm/report_array.sbatch
```

Each per-slide page has:

- **Point-based panels (always)** -- from the warped nuclei + `qc.json`: DAPI-vs-warped-H&E
  centroid scatter coloured matched/unmatched; a matched-pair displacement quiver; a histogram
  of matched nucleus offsets (um) with the `median_um` line; the two binned nucleus-density maps
  (DAPI, warped H&E) on the **same bins** `density_r` uses; and a text banner (`density_r`,
  `median_um`, `occupancy`, negative-control `density_collapse`, `status`, `chosen` + rule).
- **Raster panels (only if `registered/*.ome.tif*` exists)** -- a magenta/green DAPI-vs-H&E
  overlay and a `SimpleITK` checkerboard mosaic. Skipped cleanly if the warped image hasn't been
  generated (or for `coarse`/`rescued` slides without a raster -- those render point panels only).
- **Annotation panels (only if `cell_labels.parquet` exists)** -- cells coloured by transferred
  region label, plus a background-cell panel (cells landing off tissue = a registration smell).

It handles `status == "no_nuclei"` and all four `chosen` states (`micro` / `nomicro` / `coarse`
/ `rescued`) without crashing. `SimpleITK` is optional -- if absent, the checkerboard falls back
to a numpy mosaic.

**Human scoring (optional):** `notebooks/review_alignment.ipynb` is a thin viewer that flips
through the generated `report.png` panels and writes `output/alignment_validation.csv`
(`slide, chosen, aligned [0/1/2], has_fold, note`). It only reads the batch outputs and collects
scores -- all computation stays in `report.py`.

### Cohort QC overlay PDF (`run_cohort_qc.py`)

A single **cohort** PDF that answers "did every slide register?" at a glance: **page 1** is a
summary table of the **selected/final protocol per slide** + `density_r` / `median_um` +
disposition (accepted / manual-review / rescued, from the real `provenance.gate`); then **one page
per slide** shows the H&E nuclei (red) overlaid on the Xenium DAPI nuclei (grey) -- for the
**selected protocol only** (the losing protocol and any failed/rescued intermediate are not drawn).
Runs on results, no re-registration.

```bash
python run_cohort_qc.py --samples samples.csv --config config.json   # -> output/cohort_qc.pdf
python run_cohort_qc.py --samples samples.csv --config config.json --out mycohort.pdf
```

Reads each sample's `qc.json` decision + `he_nuclei_<chosen>.npy` + the Xenium cells; slides
without a decision/nuclei are listed and skipped, never faked. The render core
(`render_cohort_pdf`) is importable so a legacy-format cohort can build the slide list itself and
reuse the identical layout. Covered by `tests/test_run_cohort_qc.py`.

## Self-healing: coarse-alignment fallback

VALIS feature matching can lock onto a wrong solution when the H&E is grossly mis-oriented
(e.g. a 90/180/270-degree rotation or a mirror). The slide then returns a **negative density-r**:
locally nothing coincides even though the footprints roughly overlap. No rigid / non-rigid /
reflection inside VALIS fixes it, because the starting orientation is wrong.

`run_qc.py` detects this (selected density-r below `COARSE_TRIGGER`, default 0.10) and runs
`coarse_align` automatically: it searches rotation x flip and, for each, finds the best
translation by FFT phase correlation, scoring by nuclei-density agreement. If it beats the
failed registration it is selected (`rule = coarse_rescue_negative_density_r`) and saved as
`he_nuclei_coarse.npy`. (On one internal slide -- a 270-degree rotation -- this moved density-r
from -0.13 to +0.76 with no manual input; treat that as a single illustrative observation, not a
benchmark.) Disable with `--no-coarse-fallback`.

The coarse fallback fixes the nuclei, QC and annotation. To turn that into a FULL registration
(proper micro / no-micro + a warpable image), `run_rescue.py` finishes the job automatically: it
picks the cardinal rotation (`coarse_align.cardinal_rotation`), losslessly pre-rotates the H&E
(`registration.prerotate_he`), re-registers with VALIS + micro, re-QCs, and -- if it beats the
coarse density-r -- adopts it (`rule = prerotate_reregister`). Add `--warp-image` to also emit
the rescued WSI.

```bash
python run_rescue.py --samples samples.csv --config config.json --sample SLIDE_A --warp-image
```

So a grossly mis-oriented slide goes: failed VALIS (negative density-r) -> auto coarse rescue
(QC + annotations usable) -> optional `run_rescue` (full re-registration + WSI), with no manual
landmark clicking.

**Known limitation -- `run_rescue` handles cardinal ROTATIONS only (0/90/180/270).** The lossless
pre-rotation (`registration.prerotate_he`, pyvips rot90/180/270) cannot express a mirror, and the
re-register step does not flip. So a slide whose mis-orientation is (or includes) a **flip** is
still *QC-rescued* -- `coarse_align` does search flips, so its `he_nuclei_coarse` and the per-cell
annotations are correct -- but `run_rescue` will not beat the coarse density-r, so it keeps
`chosen=coarse` and does **not** emit a registered image for that slide. Flipped slides therefore
get QC + annotations but no warped WSI; handle those manually (or with a flip-capable warp).

> **Planned: flip-capable orientation rescue.** A cohort whose H&E is scanned 90-deg-rotated AND
> mirrored vs the DAPI (seen in TNBC) needs the lossless pre-orientation to include a `pyvips`
> `fliphor`, and the orientation should be picked by **nucleus-coincidence median_um** rather than
> density-r (density-r is nearly mirror-insensitive, so a wrong flip can score high while every
> cell lands mirror-imaged). When that lands it must keep the single triage schema: write the
> outcome into `qc.json` (`chosen="rescued"`, `metrics["rescued"]`, `decision.sel_*`, plus
> `decision.prerotate_flip`), not a separate metrics file -- the report and `provenance.py` then
> pick it up unchanged.

## Layout

```
hest_valis/        registration, segment, concordance, select, coarse_align, annotate, xenium, report, config
run_segment.py     step 1  (StarDist env)   H&E nuclei
run_register.py    step 2  (valis env)      register + warp nuclei, both protocols
run_qc.py          step 3  (QC env)         metrics + per-slide selection + coarse fallback
run_annotate.py            (QC env)         per-cell annotation transfer
run_rescue.py              (valis env)      pre-rotate + re-register a coarse-flagged slide
run_wsi.py                 (valis env)      warp the chosen H&E image -> OME-TIFF
run_report.py              (QC env, CPU)    per-slide + cohort alignment QC report (PDF)
run_provenance.py          (QC env, CPU)    acceptance gate + move/audit table (provenance.csv)
run_select.py      aggregate decisions -> decision table + WSI manifest
slurm/             SLURM array wrappers (qc / wsi / report)
notebooks/         review_alignment.ipynb  thin human-scoring viewer over the reports
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
- **Pixel size** -> set `pixel_um` to your DAPI um/pixel
  ([Xenium is 0.2125](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs);
  centroids are in microns, so pixels = microns / pixel size, origin top-left).

## Credits

This pipeline stands on these works; please cite them if you use it (see also `CITATION.cff`).

- **VALIS** -- the whole-slide image registration engine.
  Gatenbee et al., "Virtual alignment of pathology image series for multi-gigapixel whole slide
  images," *Nature Communications* 14, 4502 (2023).
  [doi:10.1038/s41467-023-40218-9](https://doi.org/10.1038/s41467-023-40218-9) -
  [github.com/MathOnco/valis](https://github.com/MathOnco/valis)
- **HEST / HEST-1k (Mahmood Lab)** -- the `register_dapi_he` recipe this builds on.
  Jaume et al., "HEST-1k: A Dataset for Spatial Transcriptomics and Histology Image Analysis,"
  *NeurIPS 2024* (Datasets & Benchmarks, Spotlight).
  [arXiv:2406.16192](https://arxiv.org/abs/2406.16192) -
  [github.com/mahmoodlab/HEST](https://github.com/mahmoodlab/HEST) - license **CC BY-NC-SA 4.0**.
- **StarDist** -- H&E nuclei detection.
  Schmidt et al., "Cell Detection with Star-convex Polygons," *MICCAI 2018*, pp. 265-273
  [doi:10.1007/978-3-030-00934-2_30](https://doi.org/10.1007/978-3-030-00934-2_30); and for the
  H&E model, Weigert & Schmidt, "Nuclei Instance Segmentation and Classification in Histopathology
  Images with StarDist," *ISBIC 2022*.
  [github.com/stardist/stardist](https://github.com/stardist/stardist).
  The `2D_versatile_he` model was trained on MoNuSeg 2018 + TNBC (Naylor et al. 2018), which
  bounds where it generalises -- expect to validate (or retrain) on tissue unlike that.
- **Xenium outputs / coordinate frame** -- 10x Genomics,
  [Understanding Xenium Outputs](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs)
  (origin top-left; 0.2125 um/px; the post-Xenium H&E is on a different microscope and not
  pre-registered).
- **BioFormats** ([openmicroscopy.org/bio-formats](https://www.openmicroscopy.org/bio-formats/))
  and **OME-TIFF** ([ome-model docs](https://ome-model.readthedocs.io/en/stable/ome-tiff/)).

## License

**CC BY-NC-SA 4.0** (see `LICENSE`). This pipeline is built on the HEST recipe/resource, which is
CC BY-NC-SA 4.0, so the repo adopts the same license to respect those terms. **Downstream and
commercial use is constrained by the upstream (HEST) terms: non-commercial, share-alike.** VALIS
and StarDist carry their own licenses -- consult those projects for theirs.
