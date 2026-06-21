"""Tests for hest_valis.registration.ometiff_pages -- the truncated-WSI guard.

A warp killed mid-write (preemption / walltime / OOM) leaves a multi-GB OME-TIFF whose
BigTIFF page table was never finalised: the first-IFD offset is still the 0 placeholder, so
no reader can enumerate a single page. ometiff_pages must report 0 for such a corpse and >0
for a real image -- that is what lets warp_image refuse to publish a truncated file and lets
run_wsi avoid skipping one as 'done'.

registration.py imports valis_hest (JVM/BioFormats) lazily, so importing it here needs no
registration env -- only tifffile/numpy.
"""
import os
import struct
import sys

import numpy as np
import tifffile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hest_valis.registration import ometiff_pages


def test_valid_ometiff_has_pages(tmp_path):
    p = tmp_path / "good.ome.tiff"
    tifffile.imwrite(str(p), np.zeros((16, 16), dtype="uint8"))
    assert ometiff_pages(str(p)) > 0


def test_truncated_bigtiff_reports_zero(tmp_path):
    # the exact real-world corruption: a BigTIFF header whose first-IFD offset is still the
    # 0 placeholder (writer killed before finalising the page table). Bytes on disk, no pages.
    p = tmp_path / "truncated.ome.tiff"
    header = (b"II"                       # little-endian
              + struct.pack("<H", 43)     # BigTIFF magic (0x002B)
              + struct.pack("<H", 8)      # bytesize of offsets
              + b"\x00\x00"               # constant
              + struct.pack("<Q", 0))     # first-IFD offset NEVER written -> still 0
    p.write_bytes(header + b"\x00" * 4096)
    assert ometiff_pages(str(p)) == 0


def test_missing_file_reports_zero(tmp_path):
    assert ometiff_pages(str(tmp_path / "does_not_exist.ome.tiff")) == 0
