import argparse, pathlib
import numpy as np
import os,json

import pygama.lgdo.lh5_store as lh5
import pygama.math.peak_fitting as pgf
import pygama.pargen.AoE_cal as cae
import pygama.pargen.ecal_th as ecal

import time
import logging


def build_hit(f_dsp:str, hit_dict:dict, f_hit:str,  
                copy_cols:list=None, overwrite:bool=True):

    t_start = time.time()

    if os.path.exists(f_hit) and overwrite==True:
        os.remove(f_hit)

    for channel in hit_dict:
        log.info(f'Starting Channel: {channel}')
        pars_dict = hit_dict[channel]
        log.debug(f'Using {pars_dict}')
        
        if len(pars_dict )==0:
            log.debug(f'No hit pars foudn for channel {channel}')
            pass
        else:
        
            data = lh5.load_dfs(f_dsp, ["cuspEmax", "trapEmax", "zacEmax", "dt_eff", "A_max"], f'{channel}/dsp')

            if "cusp_ctc" in pars_dict:
                data["cusp_ctc"] = ecal.apply_ctc(data["cuspEmax"],data["dt_eff"], pars_dict["cusp_ctc"]["pars"]["a"])
                data["zac_ctc"] = ecal.apply_ctc(data["zacEmax"],data["dt_eff"], pars_dict["zac_ctc"]["pars"]["a"])
                data["trap_ctc"] = ecal.apply_ctc(data["trapEmax"],data["dt_eff"], pars_dict["trap_ctc"]["pars"]["a"])

                data["cusp_cal"] = pgf.poly(data["cusp_ctc"], [pars_dict["Energy_cal_cuspEmax_ctc"]["pars"]["a"],pars_dict["Energy_cal_cuspEmax_ctc"]["pars"]["b"]])
                data["zac_cal"] = pgf.poly(data["zac_ctc"], [pars_dict["Energy_cal_zacEmax_ctc"]["pars"]["a"],pars_dict["Energy_cal_zacEmax_ctc"]["pars"]["b"]])
                data["trap_cal"] = pgf.poly(data["trap_ctc"], [pars_dict["Energy_cal_trapEmax_ctc"]["pars"]["a"],pars_dict["Energy_cal_trapEmax_ctc"]["pars"]["b"]])

                data["AoE"] = data["A_max"]/pgf.poly(data["cuspEmax"], [pars_dict["Energy_cal_cuspEmax_ctc"]["pars"]["a"],pars_dict["Energy_cal_cuspEmax_ctc"]["pars"]["b"]])
            else:
                data["cusp_cal"] = pgf.poly(data["cuspEmax"], [pars_dict["Energy_cal_cuspEmax"]["pars"]["a"],pars_dict["Energy_cal_cuspEmax"]["pars"]["b"]])
                data["zac_cal"] = pgf.poly(data["zacEmax"], [pars_dict["Energy_cal_zacEmax"]["pars"]["a"],pars_dict["Energy_cal_zacEmax"]["pars"]["b"]])
                data["trap_cal"] = pgf.poly(data["trapEmax"], [pars_dict["Energy_cal_trapEmax"]["pars"]["a"],pars_dict["Energy_cal_trapEmax"]["pars"]["b"]])

                data["AoE"] = data["A_max"]/pgf.poly(data["cuspEmax"], [pars_dict["Energy_cal_cuspEmax"]["pars"]["a"],pars_dict["Energy_cal_cuspEmax"]["pars"]["b"]])

            

            aoe_mu_pars = [pars_dict["AoE_Classifier"]["pars"]["a"],pars_dict["AoE_Classifier"]["pars"]["b"]]
            aoe_sigma_pars = [pars_dict["AoE_Classifier"]["pars"]["c"],pars_dict["AoE_Classifier"]["pars"]["d"]]

            data["AoE_classifier"] = cae.get_classifier(data["AoE"], data["cusp_cal"], aoe_mu_pars, aoe_sigma_pars)
            data['Passed_AoE_low_cut'] = (data["AoE_classifier"]>pars_dict["AoE_Low_Cut"]["pars"]["a"])
            data['Passed_AoE_double_cut'] = (data["AoE_classifier"]>pars_dict["AoE_Low_Cut"]["pars"]["a"])&(data["AoE_classifier"]<pars_dict["AoE_Double_Sided_Cut"]["pars"]["b"])

            log.info(f'Writing to file: {f_hit}')


            out_cols = ["AoE_classifier", 'Passed_AoE_low_cut', 'Passed_AoE_double_cut', "cusp_cal", "zac_cal", "trap_cal"]
            if copy_cols is not None:
                out_cols = out_cols + copy_cols
            sto=lh5.LH5Store()
            col_dict = {col : lh5.Array(data[col].values, attrs={'units':''}) for col in out_cols}
            tb_hit = lh5.Table(size=len(data), col_dict=col_dict)
            tb_name = f'{channel}/hit'
            sto.write_object(tb_hit, tb_name, f_hit)
            log.info(f'Finished Channel: {channel}')
    t_elap = time.time() - t_start
    log.info(f'Done!  Time elapsed: {t_elap:.2f} sec.')

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="input file", type=str)
argparser.add_argument("--pars_file", help="hit pars file", type=str)

argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)

argparser.add_argument("--log", help="log_file", type=str)

argparser.add_argument("--output", help="output file", type=str)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode='w')
logging.getLogger('numba').setLevel(logging.INFO)
logging.getLogger('parse').setLevel(logging.INFO)
logging.getLogger('pygama.lgdo.lh5_store').setLevel(logging.INFO)
logging.getLogger('h5py._conv').setLevel(logging.INFO)


log = logging.getLogger(__name__)

with open(args.pars_file) as f:
    pars_dict = json.load(f)


pathlib.Path(os.path.dirname(args.output)).mkdir(parents=True, exist_ok=True)
build_hit(args.input, pars_dict, f_hit =args.output, overwrite=False)
