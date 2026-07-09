from __future__ import annotations

import pytest
from dbetto import AttrsDict

from legenddataflow.scripts.flow.build_chanlist import (
    get_chanlist,
    get_par_chanlist,
    get_plt_chanlist,
)


class FakeValidity:
    """Stand-in for a pre-compiled metadata catalog with a valid_for method."""

    def __init__(self, mapping):
        self.mapping = AttrsDict(mapping)

    def valid_for(self, timestamp, system=None):  # noqa: ARG002
        return self.mapping


CHANNELMAP = FakeValidity(
    {
        "DET2": {"name": "DET2", "system": "geds"},
        "DET1": {"name": "DET1", "system": "geds"},
        "DET3": {"name": "DET3", "system": "geds"},
        "SIPM1": {"name": "SIPM1", "system": "spms"},
    }
)

DET_STATUS = FakeValidity(
    {
        "DET1": {"processable": True},
        "DET2": {"processable": True},
        "DET3": {"processable": False},
        "SIPM1": {"processable": True},
    }
)

TIMESTAMP = "20230101T123456Z"


@pytest.fixture
def setup(tmp_path):
    return {"paths": {"tmp_par": str(tmp_path / "tmp/par")}}


def test_get_chanlist():
    # sorted, restricted to the system, non-processable channels dropped
    assert get_chanlist(TIMESTAMP, "cal", DET_STATUS, CHANNELMAP, "geds") == [
        "DET1",
        "DET2",
    ]
    assert get_chanlist(TIMESTAMP, "cal", DET_STATUS, CHANNELMAP, "spms") == ["SIPM1"]


def test_get_chanlist_missing_status():
    det_status = FakeValidity({"DET1": {"processable": True}})
    with pytest.raises(RuntimeError, match="not found in the status map"):
        get_chanlist(TIMESTAMP, "cal", det_status, CHANNELMAP, "geds")


def test_get_par_chanlist(setup):
    files = get_par_chanlist(
        setup,
        f"all-l200-p03-r000-cal-{TIMESTAMP}-channels",
        "dsp",
        DET_STATUS,
        CHANNELMAP,
        "geds",
        name="eopt",
    )
    assert len(files) == 2
    assert files[0].endswith(f"l200-p03-r000-cal-{TIMESTAMP}-DET1-par_dsp_eopt.yaml")
    assert "DET2" in files[1]


def test_get_plt_chanlist(tmp_path):
    setup = {"paths": {"tmp_plt": str(tmp_path / "tmp/plt")}}
    files = get_plt_chanlist(
        setup,
        f"all-l200-p03-r000-cal-{TIMESTAMP}-channels",
        "dsp",
        DET_STATUS,
        CHANNELMAP,
        "geds",
    )
    assert len(files) == 2
    assert all(f.endswith(".pkl") for f in files)
    assert files[0].endswith(f"l200-p03-r000-cal-{TIMESTAMP}-DET1-plt_dsp.pkl")
