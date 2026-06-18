"""HEST-VALIS: align H&E whole-slide images onto Xenium DAPI, with per-slide micro selection.

Modules
  registration : register H&E onto DAPI (+optional micro), warp points / warp image
  segment      : StarDist nuclei on H&E
  concordance  : QC metrics (coincidence, density-r, occupancy, negative control)
  select       : per-slide micro vs no-micro decision rule
  xenium       : load Xenium nuclei into the DAPI frame, tissue mask
  config       : load samples.csv + config.json
"""
__all__ = ["registration", "segment", "concordance", "select", "xenium", "config"]
__version__ = "1.0"
