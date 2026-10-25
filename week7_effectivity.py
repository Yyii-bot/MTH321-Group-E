"""
Week 7 -- effectivity index for the adaptive-refinement study.

The team's adaptive loop (mark - refine) needs one piece of validation
evidence that is easy to forget: the error *estimator* has to track the
*true* error, otherwise "adaptive wins" is unfounded.  Since u* is known,
we can compute the effectivity index

    eff = eta / ||u* - u_h||_energy

on every adaptive level.  eff should stay bounded (roughly constant is
ideal); if it drifts or blows up, the estimator is not reliable and the
adaptive comparison is invalid.

This script reads the CSV that the adaptive driver should emit, one row
per adaptive level:

    level, n_dof, energy_error, eta

(energy_error can be left empty -- then this script recomputes it is NOT
possible without the mesh, so fill it in from the driver using the same
corner-aware quadrature as week6_lshape_fem.py.  Keep the quadrature
identical between the uniform and adaptive runs, or the comparison dies.)

Run:  python week7_effectivity.py path/to/adaptive_levels.csv
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def effectivity(csv_path):
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    ndof = data["n_dof"]
    err = data["energy_error"]
    eta = data["eta"]

    eff = eta / err
    print("| level | N_dof | energy error | eta | effectivity eta/err |")
    print("|---|---|---|---|---|")
    for k in range(len(ndof)):
        print(f"| {int(data['level'][k])} | {int(ndof[k])} | {err[k]:.4e} | {eta[k]:.4e} | {eff[k]:.2f} |")

    # slopes vs N_dof, for the adaptive-vs-uniform comparison figure
    print("\nobserved energy-error decay vs N_dof (pairwise):")
    for k in range(len(ndof) - 1):
        p = np.log2(err[k] / err[k + 1]) / np.log2(ndof[k + 1] / ndof[k])
        print(f"  {int(ndof[k])} -> {int(ndof[k+1])}: p = {p:.2f}   "
              f"(uniform meshes gave ~1/3; a well-designed adaptive sequence "
              f"should approach 1/2)")

    if np.any(eff > 5) or np.any(eff < 0.2):
        print("\nWARNING: effectivity index is not stable -- do not trust the "
              "estimator (and therefore not the adaptive comparison) until "
              "the indicator is fixed.")
    else:
        print("\neffectivity looks stable -- the estimator tracks the true error.")

    fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
    ax[0].loglog(ndof, err, "o-", label="true energy error")
    ax[0].loglog(ndof, eta, "s--", label="estimator eta")
    ax[0].set_xlabel("N_dof"); ax[0].set_ylabel("energy error / eta")
    ax[0].legend(); ax[0].set_title("adaptive sequence")
    ax[1].semilogx(ndof, eff, "o-")
    ax[1].set_xlabel("N_dof"); ax[1].set_ylabel("effectivity index")
    ax[1].set_title("reliability evidence")
    fig.tight_layout()
    out = "week7_effectivity.png"
    fig.savefig(out, dpi=150)
    print(f"figure saved to {out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    effectivity(sys.argv[1])
