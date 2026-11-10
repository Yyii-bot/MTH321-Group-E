# Mixed Boundary-Value Problems: MTH321 Group E

A numerical PDE project on **boundary treatment, solver accuracy and corner regularity** for the Poisson equation. This repository brings together the team's mathematical design, implementation snapshot, independent validation scripts, learning materials and final report.

The main comparison is between **one-sided** and **ghost-point** treatments of Neumann and Robin boundary conditions. Additional studies examine linear solvers, matched finite-difference/finite-element discretizations, L-shaped domains, mesh grading and pure-Neumann compatibility.

**Start here:** [Final report](MTH321-Group-E-visualization-report2/Mixed_Boundary_Value_Problems_Report_TeamE%2811%29.pdf) · [Mathematical model](MTH321-Group-E-mathematical-theory2/Mixed_BVP_Modelling_Method_Design_EN.pdf) · [Implementation snapshot](MTH321-Group-E-Algorithm-Implementation2/) · [Validation scripts](MTH321-Group-E-testing-validation2/)

## Repository guide

The five folders preserve the team's role-based deliverables. Their names are retained as supplied.

```text
.
├── README.md
├── MTH321-Group-E-Project-manager2/
│   └── Learning notes and theory/implementation reference PDFs
├── MTH321-Group-E-mathematical-theory2/
│   └── Mixed_BVP_Modelling_Method_Design_EN.pdf
├── MTH321-Group-E-Algorithm-Implementation2/
│   ├── Python modules, experiment drivers and tests
│   ├── requirements.txt
│   ├── E1–E6 CSV results and PNG figures
│   └── Run manifests and implementation audit documents
├── MTH321-Group-E-testing-validation2/
│   ├── week5_*.py through week8_*.py
│   └── Validation records, milestones and audit logs
└── MTH321-Group-E-visualization-report2/
    └── Mixed_Boundary_Value_Problems_Report_TeamE(11).pdf
```

| Folder | Purpose | Suggested starting point |
|---|---|---|
| [Project management](MTH321-Group-E-Project-manager2/) | Course learning materials and supporting notes | [Finite-difference notes](MTH321-Group-E-Project-manager2/%E4%B8%80%E7%BB%B4%E5%B7%AE%E5%88%86_%E5%85%B7%E4%BD%93%E4%BE%8B%E5%AD%90%E4%B8%8E%E4%B8%80%E8%88%AC%E5%BD%A2%E5%BC%8F.pdf), [FEM and iterative-method notes](MTH321-Group-E-Project-manager2/FEM_%E8%8A%82%E7%82%B9%E8%BF%91%E4%BC%BC_%E5%BC%BA%E5%BC%B1%E5%BD%A2%E5%BC%8F%E4%B8%8E%E8%BF%AD%E4%BB%A3%E6%96%B9%E6%B3%95pdf.pdf) |
| [Mathematical theory](MTH321-Group-E-mathematical-theory2/) | Model, boundary-row derivations, error analysis and predictions | [Modelling and method design](MTH321-Group-E-mathematical-theory2/Mixed_BVP_Modelling_Method_Design_EN.pdf) |
| [Algorithm implementation](MTH321-Group-E-Algorithm-Implementation2/) | Sparse assembly, solvers, E1–E6 experiments, tests and archived outputs | [Implementation README](MTH321-Group-E-Algorithm-Implementation2/README.md) |
| [Testing and validation](MTH321-Group-E-testing-validation2/) | Independent numerical checks and development-stage evidence | [Validation README](MTH321-Group-E-testing-validation2/README.md) |
| [Visualization and report](MTH321-Group-E-visualization-report2/) | Final project narrative, figures, tables and conclusions | [Team report](MTH321-Group-E-visualization-report2/Mixed_Boundary_Value_Problems_Report_TeamE%2811%29.pdf) |

## Mathematical problem

The core problem is

$$
-\Delta u=f\quad\text{in }\Omega=(0,1)^2.
$$

Dirichlet values are imposed on the left, bottom and top edges, including **all four corners**. The **open right edge** carries either

$$
\frac{\partial u}{\partial n}=q
\qquad\text{or}\qquad
\frac{\partial u}{\partial n}+\beta u=r.
$$

The smooth manufactured data are

$$
u^*(x,y)=e^{x+y},\qquad f(x,y)=-2e^{x+y},
$$

$$
g=u^*|_{\Gamma_D},\qquad q(y)=e^{1+y},\qquad
\beta=2,\qquad r(y)=3e^{1+y}.
$$

The same exact solution supports both boundary-condition variants. For the mixed finite-difference grid, $h=1/N$ and the free-node count is $N(N-1)$.

The collection also includes a separate pure-Dirichlet FDM/P1-FEM benchmark, L-shaped-domain studies of $r^{2/3}\sin(2\phi/3)$, and a pure-Neumann compatibility and mean-gauge test.

## Main findings in the final report

- One-sided Neumann and Robin closures give approximately **first-order global accuracy**; ghost elimination gives approximately **second-order global accuracy** on the tested smooth, aligned-grid problems.
- The centred ghost derivative is second order, while the **unscaled complete ghost boundary row** has an $O(h)$ defect. Local defect and global solution error are different quantities.
- In the reported $N=128$ ghost-Neumann tolerance sweep, $\mathrm{rtol}=10^{-10}$ is the first archived level certified to meet the **10% algebraic-to-discretization error budget**.
- SuperLU has the lowest recorded total time on the tested solver grids. The comparison includes preconditioner setup; a common residual threshold does not itself certify equal solution accuracy.
- The report's smooth P1 control gives approximately $(2,1)$ orders in continuous $(L^2,\text{energy})$ error. The corner-singular solution approaches $(4/3,2/3)$ in the observed sequence.

These statements describe the **final report**. The implementation snapshot and weekly scripts preserve earlier experimental settings; see [Reproducibility and version notes](#reproducibility-and-version-notes).

## Experiments in the implementation snapshot

This map describes the supplied **code**, whose E4 differs from the final report's P1 smooth/singular study.

| ID | Implemented experiment | Representative archived files |
|---|---|---|
| E1 | Neumann/Robin × one-sided/ghost FDM; convergence, matrix and local-row audits | `E1_core_convergence.csv`, `E1_asymptotic_fits.csv`, `E1_local_row_consistency.csv` |
| E2 | SuperLU, Jacobi, red-black Gauss–Seidel and red-black SOR; stopping, accuracy and cost checks | `E2_solver_comparison.csv`, `E2_accuracy_target_summary.csv` |
| E3 | Matched pure-Dirichlet FDM/P1-FEM on the square, with exponential and sine data | `E3_matched_FDM_FEM.csv`, `E3_stencil_equivalence.csv` |
| E4 | L-shaped **FDM**: smooth harmonic $e^x\cos y$ versus the corner-singular solution | `E4_Lshape_smooth_vs_singular.csv`, `E4_asymptotic_fits.csv` |
| E5 | Singular P1 FEM on uniform, locally red-refined and radially graded L-shaped meshes | `E5_graded_mesh_convergence.csv`, `E5_error_vs_dof_fits.csv` |
| E6 | Pure-Neumann P1 FEM: compatibility correction, constant nullspace and mean constraint | `E6_pure_Neumann_audit.csv`, `E6_small_mesh_nullspace.csv` |

The [implementation folder](MTH321-Group-E-Algorithm-Implementation2/) contains **23 CSV tables and 21 PNG figures**, stored directly in that folder in this export. Examples include the archived [Neumann convergence plot](MTH321-Group-E-Algorithm-Implementation2/E1_neumann_convergence.png) and [FDM/FEM comparison](MTH321-Group-E-Algorithm-Implementation2/E3_error_comparison.png).

## Environment

The archived implementation manifest records **Python 3.12.12** with these pinned packages:

| Package | Version |
|---|---:|
| NumPy | 2.2.5 |
| SciPy | 1.16.2 |
| pandas | 2.3.2 |
| Matplotlib | 3.10.3 |
| scikit-fem | 12.0.2 |

From the repository root, create an environment and install the supplied requirements. The activation command below is for macOS/Linux.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r MTH321-Group-E-Algorithm-Implementation2/requirements.txt
```

On Windows, activate the environment using the command appropriate to your shell before running the same `pip` installation command.

The weekly validation drivers use NumPy/SciPy; the optional effectivity plot also uses Matplotlib. The full implementation requires the complete requirements file. No TeX installation is needed to run the numerical code or read the PDFs.

## Run the independent validation scripts

These entry points support the supplied flat validation-folder layout.

```bash
cd MTH321-Group-E-testing-validation2

# Manufactured forcing, flux and Robin data
python week5_verify_analytic.py

# Analytical right-boundary-row and corner-ownership checks
python week5_boundary_rows.py
```

Additional studies can be run individually:

```bash
python week6_mixed_refinement.py
python week6_lshape_fem.py
python week8_smooth_vs_nonsmooth.py
python week8_depth_checks.py
python week8_iterative_ladder.py
```

These scripts use fixed defaults in their `__main__` blocks and mostly print tables to the terminal. They do **not** automatically update the historical Markdown evidence ledger. Solver sweeps and fine-grid iterations can take substantially longer than the initial analytical checks.

The optional effectivity driver requires a user-supplied adaptive-results CSV; one is not included:

```bash
python week7_effectivity.py path/to/adaptive_levels.csv
```

Required columns are `level`, `n_dof`, `energy_error` and `eta`. The driver writes `week7_effectivity.png` in the current directory.

## Run the main implementation

The algorithm folder is a **flattened export**. Its imports expect `mixed_bvp.experiments.*`, but the original `code/mixed_bvp/` package directories are absent. The older nested-path commands in its README therefore do not run directly from the supplied layout.

The following recipe restores a separate working copy under `_reproduction/algorithm/`, preserving the supplied source files and archived outputs.

<details>
<summary><strong>Restore the original Python package layout</strong></summary>

Run this once from the repository root. It stops if the destination already exists.

```bash
python - <<'PY'
from pathlib import Path
from shutil import copy2

source = Path("MTH321-Group-E-Algorithm-Implementation2")
target = Path("_reproduction/algorithm")
if not source.is_dir():
    raise SystemExit("Run this command from the repository root.")
if target.exists():
    raise SystemExit("_reproduction/algorithm already exists; use that copy.")

package = target / "code/mixed_bvp"
experiments = package / "experiments"
tests = target / "code/tests"
docs = target / "docs"
results = target / "results"
figures = target / "figures"
for folder in (experiments, tests, docs, results, figures):
    folder.mkdir(parents=True, exist_ok=True)

for name in ("__init__.py", "core.py", "data.py", "solvers.py", "metrics.py",
             "verification.py", "dirichlet.py", "fem.py", "pure_neumann.py"):
    copy2(source / name, package / name)
for name in ("e1_core.py", "e2_solvers.py", "e3_fdm_fem.py",
             "e4_lshape_fdm.py", "e5_graded_fem.py", "e6_pure_neumann.py",
             "common.py"):
    copy2(source / name, experiments / name)
copy2(source / "__init__ (1).py", experiments / "__init__.py")
copy2(source / "run_all.py", target / "code/run_all.py")
copy2(source / "requirements.txt", target / "requirements.txt")
for path in source.glob("test_*.py"):
    copy2(path, tests / path.name)
for name in ("Precomputation_Predictions.md", "AI_Transparency_Log.md",
             "Code_Review_and_Five_Gate_Audit.md",
             "Implementation_and_Library_Boundary.md"):
    copy2(source / name, docs / name)
for pattern, folder in (("*.csv", results), ("*.json", results),
                        ("*.png", figures)):
    for path in source.glob(pattern):
        copy2(path, folder / path.name)

print("Prepared:", target.resolve())
PY
```

The working copy has the structure expected by the driver:

```text
_reproduction/algorithm/
├── code/
│   ├── run_all.py
│   ├── mixed_bvp/
│   │   └── experiments/
│   └── tests/
├── docs/
├── results/
├── figures/
└── requirements.txt
```

</details>

After restoring the layout and installing the dependencies:

```bash
cd _reproduction/algorithm

# Inspect the available options
python code/run_all.py --help

# Smaller, selected experiments for an initial check
python code/run_all.py --quick --experiments E1 E3

# Supplied unit tests: macOS/Linux shell syntax
PYTHONPATH=code python -m unittest discover -s code/tests -v

# Full E1–E6 implementation-snapshot experiment sequence
python code/run_all.py
```

The driver also accepts `--experiments E2` or another subset of E1–E6. It generates tables in `results/`, figures in `figures/` and a fresh `results/run_manifest.json` within the working copy. Quick and selective runs overwrite the same working-copy filenames as full runs, so inspect the manifest's `quick` and `experiments` fields when interpreting outputs.

The three test modules contain **22 test methods** covering assembly, boundary rows, solver updates, geometry, FEM checks and pure-Neumann constraints. Restoring the layout enables their imports; it does not establish that a new numerical run has passed.

## Reproducibility and version notes

The collection preserves multiple development stages. Original subfolder READMEs and Markdown logs are historical documentation; some refer to files or settings outside the current export.

| Item | Supplied snapshot | Final report / interpretation |
|---|---|---|
| Core mesh sequence | Main E1 defaults: `N = 8, 16, 32, 64` | Report tables use `N = 16, 32, 64, 128, 256`; fitted values need not match exactly. |
| Solver implementation | Main implementation uses SuperLU and stationary methods; GMRES is a direct-solve cross-check | CG and GMRES+ILU comparisons are present in the validation scripts and discussed in the report, but are not implemented in the main `solvers.py`. |
| Algebraic accuracy targets | Main E2 defaults to a 1% forward-error/stagnation gate | The report's archived GMRES certificate uses a 10% algebraic/discretization budget. Residual definitions also differ between branches. |
| L-shaped study | Main E4 uses FDM and smooth harmonic `exp(x) * cos(y)` | The report's smooth/singular P1 study uses `exp(x+y)` as its smooth control. |
| Local boundary-row check | `week5_boundary_rows.py` evaluates analytical formulas at `y = 0.37` | The corrected report audits an actual common grid node at `y = 3/8`. |
| Corner localization | The weekly L-shaped script uses an older coordinate-based regional mask; some logs retain a 91–92% claim | It does not reproduce the corrected radial-region result: **64.52% of squared continuous L2 error inside `r < 0.25` at `n = 128`**. |

Further practical points:

- `run_manifest.json` and `reproduction_audit.json` in the implementation folder are **archived records** of the implementation snapshot, dated October 2026. They do not establish that the later report has been freshly reproduced from this export.
- `fem_impl_ref.py` in the validation folder uses package-relative imports and is not a standalone entry point there. Use the listed `week*.py` drivers for the independent checks.
- The theory README lists some HTML, Chinese and contribution-statement files absent from this export. The repository guide above links only files actually present.
- Timing results depend on hardware and numerical-library settings. Compare norms, mesh levels, solver stops and matrix scaling before comparing numerical tables.

## Team

Roles recorded in the project report:

| Member | Role |
|---|---|
| Yi Zhu | Project Manager |
| Zhengdao Song | Mathematical Theory |
| Xiangyu Zhang | Algorithm Implementation, Visualization and Report |
| Fangchi Shi | Algorithm Implementation |
| Rong Yin | Testing and Validation |

## Audit and transparency records

- [Implementation AI transparency log](MTH321-Group-E-Algorithm-Implementation2/AI_Transparency_Log.md)
- [Implementation code review and five-gate audit](MTH321-Group-E-Algorithm-Implementation2/Code_Review_and_Five_Gate_Audit.md)
- [Implementation/library boundary](MTH321-Group-E-Algorithm-Implementation2/Implementation_and_Library_Boundary.md)
- [Validation AI audit log](MTH321-Group-E-testing-validation2/week8_ai_audit_log.md)
- [Validation evidence ledger](MTH321-Group-E-testing-validation2/results_verified.md)

These records document development and review history. Interpret their numerical claims alongside the version notes above and the corresponding scripts or report sections.

## Licence

No repository-wide licence file is included in the supplied collection.

