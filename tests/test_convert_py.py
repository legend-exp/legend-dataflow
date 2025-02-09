import numpy as np
from legenddataflow.convert_np import convert_dict_np_to_float


def test_convert_dict_np_to_float():
    input_dict = {
        "a": np.float64(2.2),
        "b": {"c": [np.float64(5.5), 6]},
        "d": (np.float64(8.8), 9),
    }
    expected_dict = {"a": 2.2, "b": {"c": [5.5, 6]}, "d": [8.8, 9]}
    assert convert_dict_np_to_float(input_dict) == expected_dict
    test2_dict = {
        "a": np.float32(2.2),
    }
    assert np.abs(test2_dict["a"] - 2.2) < 10**-6
