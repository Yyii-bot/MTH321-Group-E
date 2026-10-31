# PDE BVP Project 2 — Mixed Boundary Data (Direction 2)

**Team E document pack** — MTH321, Project 2 (PDE BVP), Weeks 5–9

This pack contains the **modelling and method-design deliverables** for the study *Verified Finite-Difference Methods for Mixed Boundary-Value Problems* (Poisson equation on the unit square: Dirichlet data on three sides, Neumann or Robin data on the open right edge; one-sided vs ghost-point closures; direct vs iterative solvers; FDM vs P1-FEM benchmark; L-shaped-domain corner-regularity innovation).

## Contents

| File | Language | Description |
|---|---|---|
| `Mixed_BVP_Modelling_Method_Design_EN.html` / `.pdf` | English | **English edition** of the modelling & method-design document (10 sections + 3 appendices, A4, 14 pp.) |
| `混合边界BVP_建模与方法设计.html` / `.pdf` | Chinese | **Chinese original** of the same document (delivered to the team first) |
| `ICS-Project2-Song-Zhengdao.docx` / `.pdf` | English | Individual Contribution Statement + AI Transparency Log entries for Song Zhengdao (Mathematical Theory) |
| `Mixed_Boundary_Value_Problems_Report_TeamE.pdf` | English | The final team report (distributed separately via WeChat) |

## What the modelling document contains

| Section | Content |
|---|---|
| §1 Background | physical meaning of Dirichlet/Neumann/Robin data (heat exchangers, pressure-Poisson); weak-form derivation chain; why the boundary closure *is* part of the method |
| §2 Continuous model | problem statement & sign conventions; manufactured solution $u^* = e^{x+y}$ with all data; solvability (incl. pure-Neumann compatibility & gauge) |
| §3 FDM discretisation | grid & unknowns ($n = N(N-1)$); five-point stencil with exact defects; **both right-edge closures derived to matrix level** (one-sided and ghost, Neumann and Robin); matrix properties ($M$-matrix, nonsymmetry $\delta_A$); scaled residuals; convergence protocol |
| §4 Algebraic layer (Depth 2) | SuperLU vs Jacobi/Gauss–Seidel/SOR; spectral-radius theory ($\rho_J=\cos\pi h$, $\omega_{\mathrm{opt}}$); iteration-count prediction table; stopping-criterion analysis (residual vs iterate difference); E2 protocol |
| §5 FDM vs FEM benchmark (Depth 3) | pure-Dirichlet rectangular benchmark; coincident-node/DOF design; the exact stiffness identity $h^2A_{\mathrm{FDM}} = K_{P1}$; matched comparison protocol |
| §6 Regularity & degradation (Depth 4) | Céa + interpolation mechanism; $H^{1+s}$ → order-degradation rates |
| §7 L-shaped innovation | corner singularity $\lambda_1 = \pi/\omega = 2/3$; branch convention; graded-mesh scaling argument (threshold $\mu^* = 1/\lambda = 3/2$); red-refinement vs radial-grading routes |
| §8–§10 | experiment matrix E1–E6 with written predictions; 14-item risk checklist; course-protocol mapping; implementation interface (function signatures) |
| Appendices | A: boundary-row cards with $u^*$ substitution audits; B: spectral-radius derivation; C: master prediction cheat-sheet |

## Experiment matrix (summary)

| ID | Experiment | Key prediction (written before running) |
|---|---|---|
| E1 | Two closures × {Neumann, Robin} | one-sided → order 1; ghost → order 2 |
| E2 | Direct vs Jacobi/GS/SOR | step slopes −2/−2/−1; direct solver wins throughout (no crossover — report explicitly) |
| E3 | FDM vs P1-FEM, matched nodes/DOFs | both order 2; identical nonzero patterns; constant gap from load vectors |
| E4 | Smooth vs $r^{2/3}$ singular on the L-domain | smooth → 2; singular $L_\infty \to 2/3$, $L^2_h \to 4/3$ |
| E5 | Graded vs uniform meshes | uniform $N_{\mathrm{dof}}^{-1/3}$ → graded recovers $N_{\mathrm{dof}}^{-1/2}$ (energy) |
| E6 | Pure-Neumann compatibility & gauge (optional) | $A\mathbf{1} = \mathbf{0}$ before the gauge |

## Usage notes

- **HTML**: double-click to open in any browser. Formulas render via MathJax loaded from a CDN — **an internet connection is required** on first view.
- **PDF**: self-contained; formulas are fixed into the page. Use the PDF for sharing/printing (e.g. via WeChat) and the HTML as the editable source.
- **Programming stage**: §10.2 of the document specifies the implementation interface (assembly, solvers, error metrics, L-shaped and graded-mesh utilities) and a per-experiment CSV output schema.
- **Verification culture**: all predictions were written before any solve (pre-computation protocol). Keep the protocol: predictions are never edited after the run; discrepancies are reported and explained.

## Note on the English edition

The English edition incorporates **one mathematical correction** relative to the Chinese original, in §5.3 (FDM–P1 benchmark): the earlier draft predicted extra $-1$ diagonal-neighbour couplings in the P1 stiffness row (row density 7). The team's assembly audit falsified this — in a right-isosceles triangulation the two hypotenuse-end basis gradients are orthogonal, so **all** diagonal couplings vanish and $h^2A_{\mathrm{FDM}} = K_{P1}$ exactly (audit discrepancy $2.0\times10^{-15}$; report Eq (20)). The corrected identity, the resulting analysis of the error-constant gap (load moments vs nodal samples), and the audit trail are consistent with the team report and the ICS AI-transparency log (items A2/A3).

## Team

| Member | Role |
|---|---|
| Yi Zhu | Project Manager |
| Zhengdao Song | Mathematical Theory (modelling; this document pack) |
| Xiangyu Zhang | Algorithm Implementation; Visualisation and Report |
| Rong Yin | Testing and Validation |

---
*MTH321 · Department of Applied Mathematics · October 2026*


