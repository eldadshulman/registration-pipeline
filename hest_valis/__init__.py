"""HEST-VALIS: align H&E whole-slide images onto Xenium DAPI, with per-slide micro selection.

Modules
  registration : register H&E onto DAPI (+optional micro), warp points / warp image
  segment      : StarDist nuclei on H&E
  concordance  : QC metrics (coincidence, density-r, occupancy, negative control)
  select       : per-slide micro vs no-micro decision rule
  coarse_align : self-healing fallback for a failed (negative density-r) registration
  annotate     : per-cell annotation transfer (tag each Xenium cell with its H&E region)
  xenium       : load Xenium nuclei/cells into the DAPI frame, tissue mask, H&E pixel size
  report       : non-interactive per-slide + cohort alignment QC report (PDF)
  config       : load samples.csv + config.json
"""
__all__ = ["registration", "segment", "concordance", "select", "coarse_align",
           "annotate", "xenium", "report", "config"]
__version__ = "1.0"
