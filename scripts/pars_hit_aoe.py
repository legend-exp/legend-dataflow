import numpy as np
import os,json
import pathlib
import argparse

from pygama.pargen.aoe_cal import cal_aoe


argparser = argparse.ArgumentParser()
argparser.add_argument("files", help="files", nargs='*',type=str)
argparser.add_argument("--ecal_file", help="ecal_file",type=str, required=True)
argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--plot_file", help="plot_file",type=str)
argparser.add_argument("--hit_pars", help="hit_pars",type=str)
args = argparser.parse_args()

with open(args.files[0]) as f:
    files = f.read().splitlines()
    files = sorted(files)

with open(args.ecal_file, 'r') as o:
    cal_dict = json.load(o)

energy_param = 'cuspEmax'
cal_energy_param = 'cuspEmax_ctc'
    

#out_dict = cal_aoe(files, cal_dict, energy_param, cal_energy_param, dt_corr=False, cut_parameters=cut_parameters, plot_savepath=args.plot_file)

"""
cal_dict = "ecal": {
            "function": "a+b*x",
            "parameter": "cuspEmax_ctc",
            "params":{
                "a":1,
                "b":1
            }
        },
        "eres":{
            "function": "sqrt(a+b*x)",
            "parameter": "cuspEmax_ctc",
            "params":{
                "a":1,
                "b":1
            }
        },
        "aoecal": {
            "current_parameter": "A_max",
            "energy_parameter": "cuspEmax",
            "cal_energy_parameter": "cuspEmax_ctc",
            "mean":{
                "function": "a+b*x",
                "params":{
                    "a":1,
                    "b":1
                    }
                },
            "sigma":{
                "function": "sqrt(a+(B/x)**c)",
                "params":{
                    "a":1,
                    "b":1,
                    "c":1
                }
            }
        }
"""

#out_dict = {"ecal_pars":cal_dict,"aoe_pars":cal_dict}

cal_dict = {
    "001": {
        "operation1": {
            "inputs":["cuspEmax"],
            "outputs":["cal_energy"] ,
            "function":"a+b*cuspEmax",
            "pars":{"a":1, "b":1}
        },

        "operation2": {
            "inputs":["A_max", "cuspEmax", "cuspEmax_ctc"],
            "outputs":["aoe_classifier"] ,
            "function":"(((A_max/cuspEmax)/(a+b*cuspEmax_ctc)) -1 )/sqrt(c+(d/cuspEmax_ctc)**e)",
            "pars":{"a":1, "b":1, "c":1, "d":1, "e":1}
        },

        "operation3": {
            "inputs":["aoe_classifier"],
            "outputs":["aoe_pass"] ,
            "function":"a<aoe_pass<b",
            "pars":{"a":-1, "b":1}
        }
    },

    "002": {
        "operation1": {
            "inputs":["cuspEmax"],
            "outputs":["cal_energy"] ,
            "function":"a+b*cuspEmax",
            "pars":{"a":1, "b":1}
        },

        "operation2": {
            "inputs":["A_max", "cuspEmax", "cuspEmax_ctc"],
            "outputs":["aoe_classifier"] ,
            "function":"(((A_max/cuspEmax)/(a+b*cuspEmax_ctc)) -1 )/sqrt(c+(d/cuspEmax_ctc)**e)",
            "pars":{"a":1, "b":1, "c":1, "d":1, "e":1}
        },

        "operation3": {
            "inputs":["aoe_classifier"],
            "outputs":["aoe_pass"] ,
            "function":"a<aoe_pass<b",
            "pars":{"a":-1, "b":1}
        }
        
    }

}

if args.hit_pars is not None:
    pathlib.Path(os.path.dirname(args.hit_pars)).mkdir(parents=True, exist_ok=True)
    with open(args.hit_pars, 'w') as w:
        json.dump(out_dict,w, indent=4)
pathlib.Path(args.plot_file).touch()