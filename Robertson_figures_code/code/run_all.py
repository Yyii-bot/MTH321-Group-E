"""run_all.py — reproduce every experiment, table, and figure of the report.

Usage:  python run_all.py
Output: results/summary.json  (all tabulated numbers quoted in the report)
        results/*.npz         (time series behind the figures)
Then:   python make_figures.py  (produces figures/fig_*.pdf/.png)
"""
import json
import os
import time

import numpy as np

import experiments as ex
from model import output_grid

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)


def main():
    t_start = time.perf_counter()
    summary = {}

    print("[1/8] reference solution ...")
    ref = ex.reference()
    summary["reference"] = dict(tighten_diff=ref["tighten_diff"],
                                y40=[float(v) for v in ref["y40"]])
    np.savez(os.path.join(RES, "reference.npz"), t=ref["t"], y=ref["y"])

    print("[2/8] eigenanalysis ...")
    eig = ex.eigenanalysis(ref)
    summary["eigenanalysis"] = dict(max_imag=eig["max_imag"], checks=eig["checks"])
    np.savez(os.path.join(RES, "eig.npz"), t=eig["t"], lmax=eig["lmax"],
             lmin=eig["lmin"], S=eig["S"])

    print("[3/8] Prothero-Robinson order studies ...")
    summary["pr_nonstiff"] = ex.pr_order_study(-1.0)
    summary["pr_stiff"] = ex.pr_order_study(-1e4)

    print("[4/8] explicit-method threshold sweeps ...")
    thr = ex.threshold_sweeps(ref["y40"])
    summary["threshold"] = thr
    summary["blowup"] = ex.blowup_times(ref)
    for name, d in ex.threshold_trajectories().items():
        np.savez(os.path.join(RES, f"traj_{name}.npz"), t=d["t"], y=d["y"])

    print("[5/8] Robertson convergence sweeps ...")
    summary["convergence"] = ex.convergence_robertson(ref["y40"])
    summary["short_interval"] = ex.rk4_short_interval_order()

    print("[6/8] adaptive step-size study ...")
    ad = ex.adaptive_study(ref["y40"])
    summary["adaptive"] = ad["rows"]
    for key, s in ad["series"].items():
        np.savez(os.path.join(RES, f"adaptive_{key}.npz"),
                 t=s["t"], h=s["h"], err=s["err"], y=s["y"])

    print("[7/8] conservation / Newton-tolerance diagnostics ...")
    summary["newton_tol"] = ex.newton_tol_sweep(ref["y40"])
    cons = ex.conservation(ref["y40"])
    summary["conservation"] = {k: dict(ymin=v["ymin"], err40=v["err40"],
                                       max_defect=float(np.max(v["defect"])))
                               for k, v in cons.items()}
    for k, v in cons.items():
        np.savez(os.path.join(RES, f"cons_{k.replace(' ', '_')}.npz"),
                 t=v["t"], defect=v["defect"])

    print("[8/8] work-precision + y2(40) tables ...")
    summary["work_precision"] = ex.work_precision(ref["y40"])
    summary["y2_table"] = ex.y2_table(ref["y40"])

    summary["wall_total_s"] = time.perf_counter() - t_start
    with open(os.path.join(RES, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"done in {summary['wall_total_s']:.1f}s -> results/summary.json")


if __name__ == "__main__":
    main()
