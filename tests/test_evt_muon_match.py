from __future__ import annotations

import numpy as np

from legenddataflow.scripts.tier.evt import _find_matching_values_with_delay

# the value carried in the muon field_config; it happens to equal the float64
# spacing at unix timestamps of this magnitude, which is what the previous
# exact-equality implementation silently relied on
JITTER = 2.384185791015625e-07
BASE = 1.7625e9
WINDOW = int(1e9 * JITTER) * JITTER  # ~56.7 us


def _legacy(arr1, arr2, jit_delay):
    """The previous implementation, kept here as the behaviour reference."""
    matching_values = []
    for delay in np.linspace(0, int(1e9 * jit_delay), 10000) * jit_delay:
        mask = np.isin(arr1, arr2 + delay, assume_unique=True)
        matching_values.extend(arr1[mask])
    return np.unique(matching_values)


def _ticks(offsets):
    """Timestamps offset from BASE by whole float64 ticks."""
    return BASE + np.asarray(offsets, dtype=float) * np.spacing(BASE)


def test_matches_only_inside_the_window():
    muon = _ticks([0])
    # 0 and 238 ticks are inside the window, 239 is past it
    ged = _ticks([0, 1, 238, 239, 5000])

    out = _find_matching_values_with_delay(ged, muon, JITTER)

    assert np.array_equal(out, _ticks([0, 1, 238]))


def test_matching_is_one_sided():
    """A germanium trigger *before* a muon must not be flagged."""
    muon = _ticks([100])
    ged = _ticks([99, 100, 101])

    out = _find_matching_values_with_delay(ged, muon, JITTER)

    assert np.array_equal(out, _ticks([100, 101]))


def test_empty_inputs_give_no_matches():
    assert len(_find_matching_values_with_delay(np.array([]), _ticks([0]), JITTER)) == 0
    assert len(_find_matching_values_with_delay(_ticks([0]), np.array([]), JITTER)) == 0


def test_agrees_with_the_previous_implementation():
    """The rewrite must reproduce the old results exactly, not merely closely."""
    rng = np.random.default_rng(0)
    for _ in range(20):
        muon = np.unique(_ticks(rng.integers(0, 400_000, size=rng.integers(1, 6))))
        ged = np.unique(_ticks(rng.integers(0, 400_000, size=400)))
        assert np.array_equal(
            _find_matching_values_with_delay(ged, muon, JITTER),
            _legacy(ged, muon, JITTER),
        )


def test_window_scales_quadratically_with_the_jitter():
    """Documents the surprising config behaviour the refactor preserves.

    ``jitter`` is not the tolerance: the window is ``int(1e9*jitter)*jitter``.
    """
    muon = _ticks([0])
    ged = BASE + np.array([0.0, 5e-4])  # 0.5 ms after the muon

    # at the configured jitter the 0.5 ms trigger is far outside the ~57 us window
    assert len(_find_matching_values_with_delay(ged, muon, JITTER)) == 1
    # a 4x larger jitter gives a ~1 ms window, not a ~230 us one, so it matches
    assert len(_find_matching_values_with_delay(ged, muon, 1e-6)) == 2
