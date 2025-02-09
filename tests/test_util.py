from datetime import datetime
from pathlib import Path

import yaml
from legenddataflow import (
    subst_vars,
    utils,
)

testprod = Path(__file__).parent / "dummy_cycle"

with (testprod / "config.yaml").open() as r:
    setup = yaml.safe_load(r)
subst_vars(setup, var_values={"_": str(testprod)})


def test_util():
    assert utils.tier_path(setup) == str(testprod / "generated/tier")
    time = datetime.now()
    assert int(utils.unix_time(time.strftime("%Y%m%dT%H%M%SZ"))) == int(
        time.timestamp()
    )
