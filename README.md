# Team Reference File Interface

Final experiments must use data verified by the members responsible for Testing & Validation. The current `development_reference.npz` is intended only for development and validation.

The NPZ file uses pickle-free arrays:

| Key | Shape / Type | Requirements |
| --- | --- | --- |
| `t_ref` | `(N,)`, float | Strictly increasing; starts at `0` and ends at `40`; contains at least 200 positive output times. |
| `y_ref` | `(N, 3)`, float | Each row contains `[y1, y2, y3]`; all values must be finite; the first row must be `[1, 0, 0]`. |
| `metadata` | Optional scalar JSON string | Documents the source, method, tolerances, independent verification, and justification for the reported significant digits. |

The recommended time grid consists of `[0]` followed by at least 200 logarithmically spaced output points from `1e-8` to `40`. The solver's internal adaptive time grid is distinct from this output grid.

The following code only demonstrates how to export the team's existing arrays in the agreed format. It does not compute or replace the reference solution:

```python
import json
import numpy as np

# t_ref and y_ref must come from calculations verified by the team.
# Do not construct trajectories from textbook sanity-check values.
metadata = {
    "source": "team_supplied",
    "validation_note": (
        "Document the independent verification actually performed, "
        "the solver, the tolerances, and the justification for the "
        "reported significant digits here."
    ),
}

np.savez_compressed(
    "data/team_reference.npz",
    t_ref=t_ref,
    y_ref=y_ref,
    metadata=json.dumps(metadata, ensure_ascii=False),
)
```

`load_reference` only performs file, array-dimension, finiteness, time-coverage, and initial-value checks; it does not automatically certify scientific accuracy. The file's SHA-256 hash and metadata are included in the experiment manifest.

CSV/NPZ sample outputs are stored in `results/<kind>/runs/`. Each `_samples.npz` file contains:

- Common output times `t`.
- Numerical solution `Y`.
- Corresponding reference solution `reference_Y`.
- Whether the run completed the full time interval.
- The time and state at the last accepted solver-native point: `final_native_t` and `final_native_Y`.

Failed runs save only the portion covered by the solver, without extrapolation.

Each run's JSON file stores parameters and scalar statistics. Longer histories, including Newton iteration residuals and adaptive accepted/attempted steps, are stored in the corresponding `_history.npz` file with the same filename stem. The JSON field `history_file` points to this file, and `history_lengths` records the array lengths. Load these arrays with `numpy.load(..., allow_pickle=False)`. Storing long arrays separately in NPZ format avoids excessively large individual files on GitHub.

To save all solver-native points, run:

```bash
python code/run_all.py --save-trajectories
```

Full trajectories are written to `results/<kind>/trajectories/` and are not committed to Git by default.
