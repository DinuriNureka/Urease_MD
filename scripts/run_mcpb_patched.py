#!/usr/bin/env python3
"""
Run MCPB.py after patching pymsmt.mcpb.gene_pre_frcmod_file (missing pymsmtError in global scope).

Usage: same as MCPB.py, e.g.
  python3 scripts/run_mcpb_patched.py -i urease.in -s 2

Set MCPB_PY to override the MCPB.py path (default /usr/local/bin/MCPB.py).
"""
from __future__ import annotations

import os
import runpy
import sys

import pymsmt.mol.element as _el

# Urease uses residue names NI1 / NI2; pymsmt METAL_PDB only has ("NI", "NI"). Without this,
# atname "NI" is mis-read as element N (first letter), breaking ion VDW lookup (N +2).
_ni = _el.METAL_PDB[("NI", "NI")]
_el.METAL_PDB[("NI1", "NI")] = _ni
_el.METAL_PDB[("NI2", "NI")] = _ni

import pymsmt.mcpb.gene_pre_frcmod_file as _gpf
from pymsmt.exp import pymsmtError

_gpf.pymsmtError = pymsmtError

mcpb_path = os.environ.get("MCPB_PY", "/usr/local/bin/MCPB.py")
sys.argv = [mcpb_path] + sys.argv[1:]
runpy.run_path(mcpb_path, run_name="__main__")
