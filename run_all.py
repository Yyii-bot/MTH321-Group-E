#!/usr/bin/env python3
"""Reproduce every table and figure for the Algorithm Implementation deliverable.

Run from the project root:
    python code/run_all.py
    python code/run_all.py --experiments E1 E3
    python code/run_all.py --quick
"""

from __future__ import annotations
import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent
ROOT = CODE_DIR.parent
sys.path.insert(0, str(CODE_DIR))

import numpy as np
import scipy
import pandas as pd
import matplotlib
import skfem

from mixed_bvp.experiments.e1_core import run_e1
from mixed_bvp.experiments.e2_solvers import run_e2
from mixed_bvp.experiments.e3_fdm_fem import run_e3
from mixed_bvp.experiments.e4_lshape_fdm import run_e4
from mixed_bvp.experiments.e5_graded_fem import run_e5
from mixed_bvp.experiments.e6_pure_neumann import run_e6


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiments", nargs="+", default=["E1", "E2", "E3", "E4", "E5", "E6"],
                        choices=["E1", "E2", "E3", "E4", "E5", "E6"])
    parser.add_argument("--quick", action="store_true",
                        help="use smaller meshes for a fast end-to-end smoke run")
    args = parser.parse_args()

    for experiment in args.experiments:
        print(f"\n===== {experiment} =====", flush=True)
        if experiment == "E1":
            run_e1(Ns=(4, 8, 16, 32) if args.quick else (8, 16, 32, 64))
        elif experiment == "E2":
            run_e2(Ns=(8, 16, 32) if args.quick else (16, 32, 64, 128),
                   extended_N=48 if args.quick else 256,
                   tolerances=(1e-9, 1e-6) if args.quick else (1e-10, 1e-6),
                   accuracy_tolerances=((1e-9, 1e-10, 1e-11) if args.quick else
                                        (1e-10, 1e-11, 1e-12, 1e-13)))
        elif experiment == "E3":
            run_e3(Ns=(4, 8, 16, 32) if args.quick else (8, 16, 32, 64))
        elif experiment == "E4":
            run_e4(Ns=(16, 32, 64, 128) if args.quick else (32, 64, 128, 256))
        elif experiment == "E5":
            run_e5(Ns=(8, 16, 32, 64) if args.quick else (16, 32, 64, 128), kappa=6)
        elif experiment == "E6":
            run_e6(Ns=(4, 8, 16, 32))

    predictions = ROOT / "docs" / "Precomputation_Predictions.md"
    manifest = {
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "scikit-fem": skfem.__version__,
        "experiments": args.experiments,
        "quick": args.quick,
        "predictions_sha256": hashlib.sha256(predictions.read_bytes()).hexdigest(),
    }
    output = ROOT / "results" / "run_manifest.json"
    output.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"\nAll requested experiments completed. Manifest: {output}", flush=True)


if __name__ == "__main__":
    main()
