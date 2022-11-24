import numpy as np
import os,json
import pathlib
import argparse
import logging
import pickle as pkl

from util.metadata_loading import *
from util.Props import *

from pygama.pargen.ecal_th import energy_cal_th
import pygama.pargen.cuts as cts
import pygama.lgdo.lh5_store as lh5
import pygama.math.histogram as pgh

log = logging.getLogger(__name__)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--files", help="files", nargs='*',type=str)
    argparser.add_argument("--ctc_dict", help="ctc_dict", nargs="*")
    
    argparser.add_argument("--configs", help="config", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    
    argparser.add_argument("--log", help="log_file", type=str)
    
    argparser.add_argument("--plot_path", help="plot_path", type=str, required=False)
    argparser.add_argument("--save_path", help="save_path", type=str)
    argparser.add_argument("--results_path", help="results_path", type=str)
    args = argparser.parse_args()

    logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode='w')
    logging.getLogger('numba').setLevel(logging.INFO)
    logging.getLogger('parse').setLevel(logging.INFO)
    logging.getLogger('matplotlib').setLevel(logging.INFO)
    logging.getLogger('lgdo.lh5_store').setLevel(logging.INFO)

    if isinstance(args.ctc_dict, list):
        database_dic = Props.read_from(args.ctc_dict)
    else:
        with open(args.ctc_dict) as f:
            database_dic = json.load(f)

    
    hit_dict = database_dic[args.channel]["ctc_params"]

    cfg_file = os.path.join(args.configs, 'key_resolve.jsonl')
    channel_dict = config_catalog.get_config(cfg_file, args.configs, args.timestamp, args.datatype)
    channel_dict = channel_dict['snakemake_rules']['pars_hit_ecal']["inputs"]['ecal_config'][args.channel]

    with open(channel_dict,"r") as r:
        kwarg_dict = json.load(r)

    if args.plot_path:
        out_dict, result_dict,plot_dict = energy_cal_th(sorted(args.files), lh5_path=f'{args.channel}/dsp', hit_dict=hit_dict, 
                                            display=1, **kwarg_dict)  

        bins = np.arange(-500,500,1)
        bls = lh5.load_nda(sorted(args.files),["baseline"], lh5_path=f'{args.channel}/dsp')["baseline"]
        bl_array,bins = np.histogram(bls, bins=bins)
        plot_dict["baseline"] = {"bl_array":bl_array, "bins":(bins[1:]+bins[:-1])/2}
        pathlib.Path(os.path.dirname(args.plot_path)).mkdir(parents=True, exist_ok=True)
        with open(args.plot_path,"wb") as f:
            pkl.dump(plot_dict,f)
    else:
        out_dict, result_dict = energy_cal_th(sorted(args.files), lh5_path=f'{args.channel}/dsp', hit_dict=hit_dict, 
                                            **kwarg_dict)
    
    out_dict.update({"cuspEmax_cal": {
                      "expression": "a*cuspEmax+b",
                      "parameters": out_dict["cuspEmax_ctc_cal"]["parameters"]
                    }})

    with open(args.save_path,'w') as fp:
        pathlib.Path(os.path.dirname(args.save_path)).mkdir(parents=True, exist_ok=True)
        json.dump(out_dict,fp, indent=4)

    with open(args.results_path,'w') as fp:
        pathlib.Path(os.path.dirname(args.results_path)).mkdir(parents=True, exist_ok=True)
        json.dump(result_dict,fp, indent=4)