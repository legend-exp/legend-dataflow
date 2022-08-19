import numpy as np
import os,json
import pathlib
import argparse
import logging

from pygama.pargen.ecal_th import energy_cal_th
import pygama.pargen.cuts as cts
import pygama.lgdo.lh5_store as lh5
import pygama.math.histogram as pgh


def run_energy_cal(files, lh5_path = 'dsp',plot_path=None, hit_dict={}):

    energy_params = ['cuspEmax_ctc', 'zacEmax_ctc', 'trapEmax_ctc'] #

    df = lh5.load_dfs([files.replace("dsp", "raw")], ["timestamp"], lh5_path.replace("dsp", "raw"))
    df["trapTmax"] = lh5.load_nda(files, ["trapTmax"], lh5_path)["trapTmax"]
    pulser_props = cts.find_pulser_properties(df, energy="trapTmax")
    if len(pulser_props) > 0:
        out_df = cts.tag_pulsers(df, pulser_props, window=0.001)
        ids = (out_df.isPulser == 1)
        log.debug(f"pulser found: {pulser_props}")
        energy = df.trapTmax.values[~ids]
        guess_keV =  2620 / np.nanpercentile(energy, 99)
    else:
        guess_keV = None
        
    out_dict, result_dict = energy_cal_th(files, energy_params, lh5_path =lh5_path, hit_dict=hit_dict, 
                                        plot_path = plot_path, p_val=0, threshold=100, guess_keV=guess_keV)  

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

    out_dict, result_dict = run_energy_cal(sorted(args.files), lh5_path=f'{args.channel}/dsp', hit_dict=hit_dict, plot_path= args.plot_path) #cut_parameters = cut_parameters, 
    
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