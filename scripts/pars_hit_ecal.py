import numpy as np
import os,json
import pathlib
import argparse
import logging

from pygama.pargen.ecal_th import energy_cal_th


def run_energy_cal(files, lh5_path = 'dsp',plot_path=None, hit_dict={}):

    energy_params = ['cuspEmax_ctc', 'zacEmax_ctc', 'trapEmax_ctc'] #
        
    out_dict, result_dict = energy_cal_th(files, energy_params, lh5_path =lh5_path, hit_dict=hit_dict, 
                                        plot_path = plot_path, p_val=0, threshold=100)  #,cut_parameters= cut_parameters

    return out_dict, result_dict 

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("files", help="files", nargs='*',type=str)
    argparser.add_argument("--ctc_dict", help="ctc_dict", type=str)
    
    argparser.add_argument("--configs", help="config", type=str, required=True)
    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    
    argparser.add_argument("--log", help="log_file", type=str)
    
    argparser.add_argument("--plot_path", help="plot_path", type=str)
    argparser.add_argument("--save_path", help="save_path", type=str)
    argparser.add_argument("--results_path", help="results_path", type=str)
    args = argparser.parse_args()

    logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode='w')
    logging.getLogger('numba').setLevel(logging.INFO)
    logging.getLogger('parse').setLevel(logging.INFO)
    logging.getLogger('matplotlib.font_manager').setLevel(logging.INFO)
    logging.getLogger('lgdo.lh5_store').setLevel(logging.INFO)

    with open(args.ctc_dict,'r') as r:
        hit_dict = json.load(r)

    out_dict, result_dict = run_energy_cal(args.files, lh5_path=f'{args.channel}/dsp', hit_dict=hit_dict, plot_path= args.plot_path) #cut_parameters = cut_parameters, 
    
    out_dict.update({"cuspEmax_cal": {
                      "expression": "@a*cuspEmax+@b",
                      "parameters": out_dict["cuspEmax_ctc_cal"]["parameters"]
                    }})

    with open(args.save_path,'w') as fp:
        pathlib.Path(os.path.dirname(args.save_path)).mkdir(parents=True, exist_ok=True)
        json.dump(out_dict,fp, indent=4)

    with open(args.results_path,'w') as fp:
        pathlib.Path(os.path.dirname(args.results_path)).mkdir(parents=True, exist_ok=True)
        json.dump(result_dict,fp, indent=4)