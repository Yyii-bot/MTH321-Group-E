"""Output provenance and configuration handling, without rerunning long studies."""

from types import SimpleNamespace

import numpy as np
import pytest

from robertson.experiments import clean, deep_update, reference_kind


def test_loading_development_data_does_not_promote_it_to_team_reference():
    ref = SimpleNamespace(metadata={"source": "development_generated", "is_team_confirmed": False})
    assert reference_kind(ref) == "development"
    ref.metadata = {"source": "team_supplied"}
    assert reference_kind(ref) == "team"
    ref.metadata = {"status": "DEVELOPMENT ONLY"}
    assert reference_kind(ref) == "development"


def test_unknown_configuration_key_cannot_silently_change_experiment():
    with pytest.raises(ValueError, match="Unknown config key"):
        deep_update({"solver_options": {"newton_tol": 1e-12}},
                    {"solver_options": {"newton_toll": 1e-8}})


def test_missing_numeric_results_remain_null_for_portable_json():
    assert clean({"E_full": np.nan, "E_40": np.inf, "completed": np.bool_(False)}) == {
        "E_full": None, "E_40": None, "completed": False,
    }


def test_run_exports_newton_histories_as_portable_npz(tmp_path):
    import json
    from scipy.integrate import solve_ivp
    from robertson.config import default_config
    from robertson.experiments import ExperimentRunner
    from robertson.model import robertson_rhs
    from robertson.reference import Reference

    cfg = default_config()
    cfg['tf'] = 0.001
    times = np.linspace(0, 0.001, 5)
    oracle = solve_ivp(robertson_rhs, (0, 0.001), [1, 0, 0],
                       method='Radau', t_eval=times, rtol=1e-12, atol=1e-15)
    assert oracle.success
    runner = ExperimentRunner(tmp_path, cfg, Reference(times,oracle.y.T), 'development')
    row = runner.run('IE', 1e-4, 'export_check')
    record_path = tmp_path/'results/development/runs'/f"{row['run_id']}.json"
    record = json.loads(record_path.read_text(encoding='utf-8'))
    assert record['summary']['y2_40'] is None
    assert not any(isinstance(v,list) for v in record['stats'].values())
    with np.load(record_path.parent/record['history_file'],allow_pickle=False) as history:
        assert len(history['newton_iters_per_step']) == row['steps']
        assert history['newton_iters_per_step'].sum() == row['newton_iters']
        assert np.max(history['newton_residuals']) <= cfg['solver_options']['newton_tol']
