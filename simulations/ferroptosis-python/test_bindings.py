"""
Test Python bindings for ferroptosis-core.

Run: python test_bindings.py
Or:  pytest test_bindings.py -v
"""

import ferroptosis_core as fc


def test_default_params():
    params = fc.default_params()
    assert isinstance(params, dict)
    assert len(params) == 20
    assert params["death_threshold"] == 10.0
    assert params["rsl3_gpx4_inhib"] == 0.92
    assert params["fenton_rate"] == 0.02


def test_invivo_params():
    params = fc.invivo_params()
    assert params["scd_mufa_rate"] > 0
    assert params["initial_mufa_protection"] > 0.3


def test_sim_cell_returns_dict():
    result = fc.sim_cell("Persister", "RSL3", seed=42)
    assert isinstance(result, dict)
    assert "dead" in result
    assert "lp" in result
    assert "gsh" in result
    assert "gpx4" in result
    assert isinstance(result["dead"], bool)
    assert isinstance(result["lp"], float)


def test_sim_cell_determinism():
    r1 = fc.sim_cell("Persister", "RSL3", seed=42)
    r2 = fc.sim_cell("Persister", "RSL3", seed=42)
    assert r1 == r2, "Same seed should produce identical results"


def test_sim_cell_different_seeds():
    r1 = fc.sim_cell("Persister", "RSL3", seed=42)
    r2 = fc.sim_cell("Persister", "RSL3", seed=99)
    # Different seeds should generally produce different lp values
    # (not guaranteed but overwhelmingly likely)
    assert r1["lp"] != r2["lp"], "Different seeds should produce different results"


def test_sim_batch():
    stats = fc.sim_batch("Persister", "RSL3", n=1000, seed=42)
    assert isinstance(stats, dict)
    assert 0.3 < stats["death_rate"] < 0.5, f"Expected ~40%, got {stats['death_rate']}"
    assert stats["ci_low"] < stats["death_rate"] < stats["ci_high"]
    assert stats["n_dead"] > 0
    assert stats["n_cells"] == 1000
    assert stats["mean_lp"] > 0
    assert stats["mean_gsh"] > 0
    assert stats["mean_gpx4"] > 0


def test_sim_batch_param_override():
    full = fc.sim_batch("Persister", "RSL3", n=1000, seed=42)
    weak = fc.sim_batch("Persister", "RSL3", n=1000, seed=42, rsl3_gpx4_inhib=0.5)
    assert weak["death_rate"] < full["death_rate"], "Weaker drug should kill fewer cells"


def test_control_low_death_rate():
    stats = fc.sim_batch("Glycolytic", "Control", n=1000, seed=42)
    assert stats["death_rate"] < 0.05, f"Control should have very low death rate, got {stats['death_rate']}"


def test_sdt_high_kill():
    stats = fc.sim_batch("Persister", "SDT", n=1000, seed=42)
    assert stats["death_rate"] > 0.95, f"SDT should kill most persisters, got {stats['death_rate']}"


def test_bad_phenotype_raises():
    try:
        fc.sim_cell("InvalidPhenotype", "RSL3", seed=42)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "InvalidPhenotype" in str(e)


def test_bad_treatment_raises():
    try:
        fc.sim_cell("Persister", "BadTreatment", seed=42)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "BadTreatment" in str(e)


def test_bad_param_override_raises():
    try:
        fc.sim_cell("Persister", "RSL3", seed=42, nonexistent_param=1.0)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "nonexistent_param" in str(e)


if __name__ == "__main__":
    tests = [v for k, v in globals().items() if k.startswith("test_") and callable(v)]
    for test in tests:
        print(f"  {test.__name__}...", end=" ")
        test()
        print("PASS")
    print(f"\n{len(tests)} tests passed")
