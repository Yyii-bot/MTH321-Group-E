# Week 5 -- milestone notes (Testing & Validation)

## What I did this week

Ran every analytic check *before* any solver existed, so that when numbers
go wrong later we know whether to blame the maths or the code.

- `week5_verify_analytic.py` -- all green:
  - Direction 2: f = -2 exp(x+y) confirmed against a fine 5-point stencil
    (not by hand-waving derivatives); boundary data sampled on all four
    edges; q = exp(1+y) > 0 confirms the outward normal on x=1 points
    the right way; Robin data r = 3 exp(1+y) satisfies r = q + 2u*.
  - L-shape: branch convention phi = theta - pi/2 in [0, 3pi/2] unit-tested
    at five points including the nasty quadrant IV; u* = 0 on both rays
    (machine zero); outer boundary data nonzero (so the problem isn't
    trivially u = 0); Delta u* = 0 verified by stencil at points with
    r >= 0.3; u* >= 0 over 1500 random domain samples.
- `week5_boundary_rows.py` -- the four right-edge rows (first-order and
  ghost, Neumann and Robin) verified by substituting sampled u* values:
  - first-order rows: residual shrinks like h (ratio 1.99 per halving),
    as a first-order closure must.
  - ghost rows: residual shrinks like h *locally* -- this is expected
    (LeVeque FDML), the O(h^2) payoff must show up in the solution error,
    which is exactly what the Week 6 tables will test.
  - corner ownership helper: asserts the four Dirichlet corners are never
    overwritten by a boundary row. Runs against the real tag grid.

## What I learned / flag for the team

The problem pack prints the PDE as "Delta u = f" but then gives
f = -2 exp(x+y), which is consistent with -Delta u = f (the test bank's
convention). We adopt -Delta u = f. Worth a line in the report so the
clinic doesn't trip on it.

## Plan for Week 6

Four-mesh refinement tables for all four closure/BC combinations, with
boundary-vs-interior errors split out, scaled residuals, and the matrix
asymmetry numbers. Same for the L-shape on uniform meshes.

## Blocker

None.
