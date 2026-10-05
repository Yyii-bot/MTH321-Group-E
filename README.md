MTH321 Group E: Numerical Solution of the Robertson Problem

A numerical study of the Robertson chemical kinetics problem for **MTH321, Project 1 (ODE IVP)**. This repository combines mathematical derivations, solver code, independent validation, numerical results, and the final report and presentation.

The project compares explicit Euler, classical fourth-order Runge-Kutta, and implicit Euler, with particular attention to stiffness, stability, Newton iteration, adaptive step sizes, and the cost of achieving a specified accuracy.

[Final report](MTH321-Group-E-visualization-report/Robertson_ODE_Report_TeamE.pdf) · [Presentation slides](MTH321-Group-E-visualization-report/Robertson_ODE_SIide_TeamE.pdf) · [Mathematical derivation](MTH321-Group-E-mathematical-theory/Robertson_Mathematical_Derivation_EN.pdf) · [Validation report](MTH321-Group-E-testing-validation/report/testing-validation-section.md)

## Repository Guide

| Directory | Contents |
| --- | --- |
| [`MTH321-Group-E-mathematical-theory`](MTH321-Group-E-mathematical-theory/) | English mathematical derivation covering the model, conservation law, Jacobians, stiffness, and numerical methods. |
| [`MTH321-Group-E-Algorithm-Implementation`](MTH321-Group-E-Algorithm-Implementation/) | Experiment entry-point wrapper, reference-data interface specification, and tests for the intended solver API. See the completeness note below. |
| [`MTH321-Group-E-testing-validation`](MTH321-Group-E-testing-validation/) | Independent validation code, seven saved validation figures, execution logs, review notes, and milestone reports. |
| [`MTH321-Group-E-visualization-report`](MTH321-Group-E-visualization-report/) | Final report, presentation slides, and the report reproduction package, including Python sources, saved experiment data, and ten report figures. |
| [`MTH321-Group-E-figures`](MTH321-Group-E-figures/) | A separate collection of PNG/PDF figures covering concentrations, convergence, stiffness, stability, Newton tolerances, adaptivity, cost, and diagnostics. |

**Reproduction entry point:** use [`MTH321-Group-E-visualization-report/Robertson_figures_code`](MTH321-Group-E-visualization-report/Robertson_figures_code/), which contains the model, solvers, experiment drivers, requirements file, and plotting scripts.

**Completeness note:** the supplied `Algorithm-Implementation` directory does not include the `robertson` package imported by its wrapper and tests. That component cannot run independently as supplied. The report package and the testing-validation package each contain their own solver implementations; their test results and output formats should be interpreted separately.

## Problem Formulation

The Robertson system is

$$
\begin{aligned}
\frac{dy_1}{dt} &= -0.04y_1 + 10^4y_2y_3, \\
\frac{dy_2}{dt} &= 0.04y_1 - 10^4y_2y_3 - 3\times10^7y_2^2, \\
\frac{dy_3}{dt} &= 3\times10^7y_2^2,
\end{aligned}
$$

with initial condition $\mathbf{y}(0)=(1,0,0)^T$ and integration interval $0\leq t\leq40$.

The total concentration satisfies

$$
y_1(t)+y_2(t)+y_3(t)=1.
$$

The widely separated decay time scales make the system stiff. Eliminating $y_3=1-y_1-y_2$ gives a two-dimensional system for stiffness analysis. The structural zero eigenvalue associated with conservation is excluded from the stiffness ratio; the ratio is evaluated where both relevant decay modes are nonzero.

## Numerical Methods and Experiments

| Method | Role in the study |
| --- | --- |
| Explicit Euler (EE) | First-order explicit baseline and investigation of stability-limited step sizes. |
| Classical RK4 | Fourth-order explicit method, used to study accuracy and stability restrictions. |
| Implicit Euler (IE) | First-order implicit method with nonlinear equations solved by Newton iteration. |
| Adaptive implicit Euler | Step-doubling error estimation and step-size control. |
| SciPy Radau and BDF | Reference calculations and production-solver comparisons. The independent validation package also uses LSODA for reference cross-checks. |

The experiments examine:

- Full and reduced Jacobians, eigenvalue histories, and stiffness ratios.
- Stability regions, explicit-method step-size thresholds, and failure behavior.
- Convergence on the Robertson and Prothero-Robinson problems.
- Newton tolerance sensitivity, conservation defects, and non-negativity.
- Adaptive step sizes, accepted/rejected steps, and tolerance sensitivity.
- Work-precision curves and runtime comparisons at matched endpoint errors.

## Reproduce the Report Figures

Use **Python 3.10 or newer**. The report package declares `numpy>=1.26`, `scipy>=1.11`, and `matplotlib>=3.8` in its requirements file.

Starting from the repository root:

```bash
cd MTH321-Group-E-visualization-report/Robertson_figures_code
python -m pip install -r code/requirements.txt

# Rebuild Figures 1 and 3-9 from the supplied experiment data.
python code/make_figures.py

# Rebuild Figure 2 and Figure 10.
python code/make_figures_supplementary.py
```

To recompute the experiment data and then regenerate all figures, run the following from the same package directory:

```bash
python code/run_all.py
python code/make_figures.py
python code/make_figures_supplementary.py
```

`run_all.py` writes experiment arrays and `summary.json` to the package's `results/` directory. Plotting scripts write to its `figures/` directory. Figures 1-9 are supplied as PNG and PDF; Figure 10, `fig_matchedcost`, is PNG only.

The supplementary script recomputes Figure 2 using a separate 100-fold tolerance refinement, from `(rtol, atol) = (1e-10, 1e-12)` to `(1e-12, 1e-14)`. For Figure 10, it uses stored benchmark timings by default. To measure timings on the current machine instead:

```bash
python code/make_figures_supplementary.py --rerun
```

Runtime measurements depend on the machine, software environment, and benchmark protocol.

## Validation

From the report package directory, run its checks with:

```bash
python code/validation_tests.py
```

These checks cover the conservation identity, analytic Jacobian, Newton solver, convergence orders, numerical conservation, Radau agreement, stiffness diagnostics, and adaptive tolerance refinement.

The separate testing-validation package has its own implementation and experiment pipeline. With the same scientific Python dependencies installed, run it from the **repository root**:

```bash
python MTH321-Group-E-testing-validation/code/validation/run_all.py
```

Its generated figures are written to `MTH321-Group-E-testing-validation/figures/validation/`, and results are printed to the console. The supplied [`run_all_output.txt`](MTH321-Group-E-testing-validation/figures/validation/run_all_output.txt) records **8/8 validation tests passing and 5/5 pipeline stages completing** for that saved run. This is evidence for the independent validation implementation, not a claim that every test in the repository has passed. Inspect the printed PASS/FAIL results when rerunning the pipeline.

## Reference Data and Output Formats

The report package constructs a Radau reference solution using `rtol=1e-12` and `atol=1e-14`, then compares it with a tighter run using `rtol=1e-13` and `atol=1e-15`. The output grid contains `t=0` plus 241 logarithmically spaced times from `1e-8` to `40`.

The saved report data use the following layouts:

| File | Contents |
| --- | --- |
| `results/reference.npz` | `t` with shape `(242,)` and `y` with shape `(3, 242)`. |
| `results/eig.npz` | Output times and eigenvalue-magnitude/stiffness histories: `t`, `lmax`, `lmin`, and `S`. |
| `results/adaptive_*.npz` | Accepted times `t`, states `y`, step sizes `h`, and error estimates `err`. |
| `results/cons_*.npz` | Times `t` and concentration-sum defects `defect`. |
| `results/summary.json` | Reference values, experiment summaries, errors, diagnostics, and timing data. |

These paths are relative to `MTH321-Group-E-visualization-report/Robertson_figures_code/`.

The [Algorithm-Implementation README](MTH321-Group-E-Algorithm-Implementation/README.md) separately specifies a team interchange format with `t_ref`, `y_ref` of shape `(N, 3)`, and optional JSON metadata. It also describes experiment manifests and sample/history exports. This specification is distinct from the saved report-package format above; files should not be treated as interchangeable without converting their keys and array orientation.

Reference refinement and cross-solver agreement provide consistency checks, not rigorous error bounds. Structural file checks alone do not establish scientific accuracy.

## Main Findings

The supplied [report](MTH321-Group-E-visualization-report/Robertson_ODE_Report_TeamE.pdf), [saved experiment summary](MTH321-Group-E-visualization-report/Robertson_figures_code/results/summary.json), and validation materials support the following observations:

- The reference endpoint is approximately $\mathbf{y}(40)=(0.7158270687,\;9.185534765\times10^{-6},\;0.2841637457)^T$.
- The stiffness ratio near $t=40$ is approximately $1.58\times10^5$.
- Frozen-Jacobian analysis suggests explicit step-size limits of about $5.9\times10^{-4}$ for EE and $8.2\times10^{-4}$ for RK4. These are local stability estimates, not guarantees for the nonlinear problem.
- Implicit Euler can remain stable at larger step sizes, but its first-order accuracy and Newton cost still matter. In the report's four matched-error timing comparisons, IE took about **5.16-5.25 times as long as EE** for the tested implementation and benchmark environment.
- Conservation, non-negativity, boundedness, and accuracy are separate diagnostics. Preserving the concentration sum alone does not demonstrate an accurate solution.

## Team

| Member | Responsibility |
| --- | --- |
| Yi Zhu | Project Manager |
| Zhengdao Song | Mathematical Theory |
| Fangchi Shi | Algorithm Implementation |
| Xiangyu Zhang | Visualization and Report |
| Rong Yin | Testing and Validation |

For a guided reading path, start with the mathematical derivation, read the final report, and then consult the reproduction package and validation materials for the supporting code and numerical evidence.
