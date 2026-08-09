from __future__ import annotations

import pytest

from legenddataflow.scripts.tier.evt import _resolve_buffer_len

PYGAMA_DEFAULT = 10**4


def test_cli_value_wins_over_the_rule_config():
    assert _resolve_buffer_len(500, {"options": {"buffer_len": 20}}) == 500


def test_falls_back_to_the_rule_config_then_the_default():
    assert _resolve_buffer_len(None, {"options": {"buffer_len": 20}}) == 20
    assert _resolve_buffer_len(None, {"options": {}}) == PYGAMA_DEFAULT
    assert _resolve_buffer_len(None, {}) == PYGAMA_DEFAULT


@pytest.mark.parametrize("bad", [0, -1, -(10**6)])
def test_non_positive_chunk_sizes_are_rejected(bad):
    """A chunk size of 0 or less is meaningless and fails obscurely later."""
    with pytest.raises(ValueError, match="must be positive"):
        _resolve_buffer_len(bad, {})
    with pytest.raises(ValueError, match="must be positive"):
        _resolve_buffer_len(None, {"options": {"buffer_len": bad}})


def test_the_message_names_where_the_value_came_from():
    with pytest.raises(ValueError, match=r"--buffer-len"):
        _resolve_buffer_len(0, {})
    with pytest.raises(ValueError, match="tier_evt rule config"):
        _resolve_buffer_len(None, {"options": {"buffer_len": 0}})


@pytest.mark.parametrize("bad", ["abc", None, [1]])
def test_non_integer_config_values_are_rejected(bad):
    if bad is None:  # `buffer_len: null` means "unset", so the default applies
        assert (
            _resolve_buffer_len(None, {"options": {"buffer_len": None}})
            == PYGAMA_DEFAULT
        )
        return
    with pytest.raises(ValueError, match="is not an integer"):
        _resolve_buffer_len(None, {"options": {"buffer_len": bad}})
