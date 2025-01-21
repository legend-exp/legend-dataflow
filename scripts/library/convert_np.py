import numpy as np


def convert_dict_np_to_float(dic):
    for key in dic:
        if isinstance(dic[key], dict):
            convert_dict_np_to_float(dic[key])
        elif isinstance(dic[key], (np.float32, np.float64)):
            dic[key] = float(dic[key])
        elif isinstance(dic[key], (list, tuple)):
            dic[key] = [
                float(x) if isinstance(x, (np.float32, np.float64)) else x for x in dic[key]
            ]
    return dic
