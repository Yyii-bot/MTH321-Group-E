"""Reference provenance, interchange validation and interpolation contracts."""

import json

import numpy as np
import pytest
from scipy.interpolate import CubicHermiteSpline

from robertson.diagnostics import align_to_reference, error_diagnostics
from robertson.model import robertson_rhs
from robertson.reference import Reference, create_development_reference, load_reference, save_reference


def _structurally_valid_data():
    # Synthetic fixture only for the file-format validator, not an ODE oracle.
    t = np.r_[0.0, np.geomspace(1e-8, 40.0, 200)]
    Y = np.zeros((len(t), 3))
    Y[:, 0] = 1.0
    return t, Y


@pytest.mark.parametrize("problem", ["missing_key", "wrong_shape", "nonfinite", "unsorted", "short_grid", "wrong_coverage", "wrong_initial", "metadata_list", "metadata_invalid_json"])
def test_reject_malformed_reference_files(tmp_path, problem):
    t, Y = _structurally_valid_data()
    data = {"t_ref": t, "y_ref": Y}
    if problem == "missing_key":
        del data["y_ref"]
    elif problem == "wrong_shape":
        data["y_ref"] = Y.T
    elif problem == "nonfinite":
        Y[5, 1] = np.nan
    elif problem == "unsorted":
        t[5] = t[4]
    elif problem == "short_grid":
        data["t_ref"], data["y_ref"] = t[::2], Y[::2]
    elif problem == "wrong_coverage":
        t[-1] = 39.0
    elif problem == "wrong_initial":
        Y[0] = [0.9, 0.1, 0.0]
    elif problem == "metadata_list":
        data["metadata"] = json.dumps(["not", "a", "mapping"])
    elif problem == "metadata_invalid_json":
        data["metadata"] = "{malformed"
    path = tmp_path / "malformed.npz"
    np.savez(path, **data)
    with pytest.raises(ValueError):
        load_reference(path)


def test_development_provenance_survives_save_and_load(tmp_path):
    original = create_development_reference(n_points=200)
    loaded = load_reference(save_reference(original, tmp_path / "development.npz"))
    np.testing.assert_array_equal(loaded.t, original.t)
    np.testing.assert_array_equal(loaded.Y, original.Y)
    assert loaded.metadata["source"] == "development_generated"
    assert loaded.metadata["is_team_confirmed"] is False
    assert "DEVELOPMENT ONLY" in loaded.metadata["status"]
    assert set(loaded.metadata["solver_stats"]) == {"Radau_base", "Radau_tight", "BDF_base", "BDF_tight"}
    assert loaded.metadata["effective_rtol_tightening"] == pytest.approx(100)
    # Cross-method agreement is a check, not a certified oracle-error bound.
    assert loaded.metadata["validation"]["tight_Radau_BDF_max_diff"] < 1e-9
    with pytest.raises(ValueError, match="extrapolation"):
        loaded.evaluate(40.01)


def test_selected_interval_hermite_equals_full_trajectory_hermite():
    t = np.linspace(0, 0.001, 101)
    Y = np.column_stack((1 - 0.04 * t, 0.04 * t, np.zeros_like(t)))
    query = np.array([0.0, 0.0001223, 0.0008841, 0.001])
    actual, rhs_calls = align_to_reference(t, Y, query)
    full = CubicHermiteSpline(t, Y, [robertson_rhs(ti, yi) for ti, yi in zip(t, Y)])
    # Interior interpolation is identical; stored nodes deliberately retain
    # their exact values instead of cubic-evaluation roundoff.
    np.testing.assert_array_equal(actual[1:3], full(query)[1:3])
    np.testing.assert_array_equal(actual[[0, 3]], Y[[0, -1]])
    assert rhs_calls <= 2 * len(query)


def test_failed_prefix_has_only_partial_error():
    ref = Reference(
        np.array([0.0, 20.0, 40.0]),
        np.array([[1.0, 0, 0], [0.9, 0, 0.1], [0.8, 0, 0.2]]),
    )
    diagnostics = error_diagnostics(ref.t[:2], ref.Y[:2], ref, completed=False)
    assert np.isnan(diagnostics["E_full"])
    assert np.isnan(diagnostics["E_40"])
    assert diagnostics["E_full_partial"] == 0.0
    assert diagnostics["comparison_points"] == 2
    assert not diagnostics["accurate_flag"]
    assert not diagnostics["completed"]


def test_accuracy_is_independent_of_conservation_and_positivity():
    ref = Reference(
        np.array([0.0, 20.0, 40.0]),
        np.array([[1.0, 0, 0], [0.9, 0, 0.1], [0.8, 0, 0.2]]),
    )
    incorrect_but_positive = np.tile([1.0, 0.0, 0.0], (3, 1))
    diagnostics = error_diagnostics(ref.t, incorrect_but_positive, ref, accuracy_threshold=1e-4)
    assert diagnostics["E_cons"] == 0.0
    assert diagnostics["nonnegative_flag"]
    assert not diagnostics["accurate_flag"]
    assert diagnostics["completed"]
