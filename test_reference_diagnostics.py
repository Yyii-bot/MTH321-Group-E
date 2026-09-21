"""No missing endpoint may be disguised as a completed experiment."""

import numpy as np
import pytest

from robertson.diagnostics import error_diagnostics
from robertson.reference import Reference, create_development_reference, load_reference, save_reference
from robertson.stability import eigenanalysis


@pytest.fixture(scope="module")
def development_reference():
    return create_development_reference(n_points=200)


def test_development_reference_records_independent_validation(development_reference):
    ref = development_reference
    assert ref.t[0] == 0.0
    assert ref.t[-1] == 40.0
    assert np.count_nonzero(ref.t > 0) >= 200
    assert ref.metadata["source"] == "development_generated"
    assert not ref.metadata["is_team_confirmed"]
    assert set(ref.metadata["solver_stats"]) == {
        "Radau_base", "Radau_tight", "BDF_base", "BDF_tight"
    }
    assert ref.metadata["validation"]["tight_Radau_BDF_max_diff"] < 1e-8
    assert np.max(np.abs(ref.Y.sum(axis=1) - 1)) < 1e-10


def test_reference_roundtrip_preserves_samples_and_provenance(development_reference, tmp_path):
    path = save_reference(development_reference, tmp_path / "reference.npz")
    restored = load_reference(path)
    np.testing.assert_array_equal(restored.t, development_reference.t)
    np.testing.assert_array_equal(restored.Y, development_reference.Y)
    assert restored.metadata["source"] == "development_generated"
    assert not restored.metadata["is_team_confirmed"]
    np.testing.assert_allclose(restored.evaluate(restored.t), restored.Y, atol=1e-15)
    with pytest.raises(ValueError, match="extrapolation"):
        restored.evaluate(40.01)


def test_reference_loader_rejects_missing_arrays_and_insufficient_grid(tmp_path):
    missing = tmp_path / "missing.npz"
    np.savez(missing, t_ref=[0.0, 40.0])
    with pytest.raises(ValueError, match="t_ref and y_ref"):
        load_reference(missing)
    short = tmp_path / "short.npz"
    np.savez(short, t_ref=[0.0, 40.0], y_ref=[[1, 0, 0], [0.8, 0, 0.2]])
    with pytest.raises(ValueError, match="200"):
        load_reference(short)


def test_completed_reference_samples_have_zero_comparison_error(development_reference):
    ref = development_reference
    diag = error_diagnostics(ref.t, ref.Y, ref)
    assert diag["completed"]
    assert diag["accurate_flag"]
    assert diag["nonnegative_flag"]
    assert diag["E_full"] < 1e-14
    assert diag["E_40"] < 1e-14


def test_incomplete_run_has_no_full_interval_or_endpoint_error(development_reference):
    ref = development_reference
    diag = error_diagnostics(ref.t[:60], ref.Y[:60], ref, completed=False)
    assert not diag["completed"]
    assert not diag["accurate_flag"]
    assert np.isnan(diag["E_full"])
    assert np.isnan(diag["E_40"])
    assert np.isfinite(diag["E_full_partial"])
    assert diag["comparison_points"] == 60


def test_failure_flag_overrides_apparent_endpoint_coverage(development_reference):
    ref = development_reference
    diag = error_diagnostics(ref.t, ref.Y, ref, completed=False)
    assert not diag["completed"]
    assert not diag["accurate_flag"]
    assert np.isnan(diag["E_40"])


def test_short_reference_does_not_rename_its_endpoint_E40():
    t = np.array([0.0, 0.001])
    Y = np.array([[1.0, 0.0, 0.0], [0.99996, 3e-5, 1e-5]])
    reference = Reference(t, Y)
    diag = error_diagnostics(t, Y, reference)
    assert diag["completed"]
    assert np.isnan(diag["E_40"])


def test_conservation_positivity_and_accuracy_are_independent(development_reference):
    ref = development_reference
    Y = ref.Y.copy()
    # Keep the invariant exactly while forcing one unphysical component.
    Y[-1, 1] -= 0.01
    Y[-1, 2] += 0.01
    diag = error_diagnostics(ref.t, Y, ref, accuracy_threshold=1e-5)
    assert diag["E_cons"] < 1e-10
    assert not diag["nonnegative_flag"]
    assert not diag["accurate_flag"]
    assert diag["first_negative_time"] == 40.0
    assert diag["E_40"] > 0.01


def test_nonfinite_trajectory_never_passes_accuracy(development_reference):
    ref = development_reference
    Y = ref.Y.copy()
    Y[-1, 0] = np.nan
    diag = error_diagnostics(ref.t, Y, ref)
    assert not diag["completed"]
    assert not diag["accurate_flag"]
    assert not diag["nonnegative_flag"]
    assert np.isnan(diag["E_40"])
    assert np.isinf(diag["E_cons"])


def test_eigenanalysis_excludes_zero_and_matches_reduced_system(development_reference):
    result = eigenanalysis(development_reference)
    assert np.all(result["t"] > 0)
    assert np.all(result["valid"])
    assert np.all(result["zero_count"] == 1)
    assert np.all(np.abs(result["zero_eigenvalue"]) < 1e-10)
    assert np.all(result["reduced_agrees"])
    assert np.all(np.isfinite(result["stiffness_ratio"]))
    assert np.all(result["stiffness_ratio"] > 1)
    for mode in ("slow", "fast"):
        np.testing.assert_allclose(
            result[f"lambda_{mode}"], result[f"reduced_lambda_{mode}"],
            rtol=1e-7, atol=1e-8,
        )


def test_bad_zero_threshold_does_not_fabricate_stiffness_ratio(development_reference):
    result = eigenanalysis(development_reference, zero_threshold=1e6)
    assert not np.any(result["valid"])
    assert np.all(np.isnan(result["stiffness_ratio"]))
