import argparse
import logging
import pickle as pkl
import time
import warnings
from pathlib import Path

import lgdo.lh5 as lh5
import numpy as np
import pygama.pargen.energy_optimisation as om  # noqa: F401
import sklearn.gaussian_process.kernels as ker
from dspeed.units import unit_registry as ureg
from legendmeta import LegendMetadata
from legendmeta.catalog import Props
from pygama.math.distributions import hpge_peak
from pygama.pargen.dsp_optimize import (
    BayesianOptimizer,
    run_bayesian_optimisation,
    run_one_dsp,
)

warnings.filterwarnings(action="ignore", category=RuntimeWarning)
warnings.filterwarnings(action="ignore", category=np.RankWarning)


argparser = argparse.ArgumentParser()

argparser.add_argument("--peak_file", help="tcm_filelist", type=str, required=True)
argparser.add_argument("--decay_const", help="decay_const", type=str, required=True)
argparser.add_argument("--inplots", help="in_plot_path", type=str)

argparser.add_argument("--log", help="log_file", type=str)
argparser.add_argument("--configs", help="configs", type=str, required=True)
argparser.add_argument("--metadata", help="metadata", type=str, required=True)

argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
argparser.add_argument("--channel", help="Channel", type=str, required=True)

argparser.add_argument("--final_dsp_pars", help="final_dsp_pars", type=str, required=True)
argparser.add_argument("--qbb_grid_path", help="qbb_grid_path", type=str)
argparser.add_argument("--plot_path", help="plot_path", type=str)

argparser.add_argument("--plot_save_path", help="plot_save_path", type=str, required=False)
args = argparser.parse_args()

logging.basicConfig(level=logging.DEBUG, filename=args.log, filemode="w")
logging.getLogger("numba").setLevel(logging.INFO)
logging.getLogger("parse").setLevel(logging.INFO)
logging.getLogger("lgdo").setLevel(logging.INFO)
logging.getLogger("h5py").setLevel(logging.INFO)
logging.getLogger("matplotlib").setLevel(logging.INFO)
logging.getLogger("dspeed.processing_chain").setLevel(logging.INFO)
logging.getLogger("legendmeta").setLevel(logging.INFO)


log = logging.getLogger(__name__)
sto = lh5.LH5Store()
t0 = time.time()

meta = LegendMetadata(path=args.metadata)
channel_dict = meta.channelmap(args.timestamp, system=args.datatype)
channel = f"ch{channel_dict[args.channel].daq.rawid:07}"

conf = LegendMetadata(path=args.configs)
configs = conf.on(args.timestamp, system=args.datatype)
dsp_config = configs["snakemake_rules"]["pars_dsp_eopt"]["inputs"]["processing_chain"][
    args.channel
]
opt_json = configs["snakemake_rules"]["pars_dsp_eopt"]["inputs"]["optimiser_config"][args.channel]

opt_dict = Props.read_from(opt_json)
db_dict = Props.read_from(args.decay_const)

if opt_dict.pop("run_eopt") is True:
    peaks_kev = np.array(opt_dict["peaks"])
    kev_widths = [tuple(kev_width) for kev_width in opt_dict["kev_widths"]]

    kwarg_dicts_cusp = []
    kwarg_dicts_trap = []
    kwarg_dicts_zac = []
    for peak in peaks_kev:
        peak_idx = np.where(peaks_kev == peak)[0][0]
        kev_width = kev_widths[peak_idx]

        kwarg_dicts_cusp.append(
            {
                "parameter": "cuspEmax",
                "func": hpge_peak,
                "peak": peak,
                "kev_width": kev_width,
                "bin_width": 5,
            }
        )
        kwarg_dicts_zac.append(
            {
                "parameter": "zacEmax",
                "func": hpge_peak,
                "peak": peak,
                "kev_width": kev_width,
                "bin_width": 5,
            }
        )
        kwarg_dicts_trap.append(
            {
                "parameter": "trapEmax",
                "func": hpge_peak,
                "peak": peak,
                "kev_width": kev_width,
                "bin_width": 5,
            }
        )

    peaks_rounded = [int(peak) for peak in peaks_kev]
    peaks = sto.read(f"{channel}/raw", args.peak_file, field_mask=["peak"])[0]["peak"].nda
    ids = np.isin(peaks, peaks_rounded)
    peaks = peaks[ids]
    idx_list = [np.where(peaks == peak)[0] for peak in peaks_rounded]

    tb_data = sto.read(f"{channel}/raw", args.peak_file, idx=ids)[0]

    t1 = time.time()
    log.info(f"Data Loaded in {(t1-t0)/60} minutes")

    if isinstance(dsp_config, (str, list)):
        dsp_config = Props.read_from(dsp_config)

    dsp_config["outputs"] = ["tp_99", "tp_0_est", "dt_eff"]

    init_data = run_one_dsp(tb_data, dsp_config, db_dict=db_dict, verbosity=0)
    full_dt = (init_data["tp_99"].nda - init_data["tp_0_est"].nda)[idx_list[-1]]
    flat_val = np.ceil(1.1 * np.nanpercentile(full_dt, 99) / 100) / 10

    if flat_val < 1.0:
        flat_val = 1.0
    elif flat_val > 4:
        flat_val = 4
    flat_val = f"{flat_val}*us"

    db_dict["cusp"] = {"flat": flat_val}
    db_dict["zac"] = {"flat": flat_val}
    db_dict["etrap"] = {"flat": flat_val}

    tb_data.add_column("dt_eff", init_data["dt_eff"])

    dsp_config["processors"].pop("dt_eff")

    dsp_config["outputs"] = ["zacEmax", "cuspEmax", "trapEmax", "dt_eff"]

    kwarg_dict = [
        {
            "peak_dicts": kwarg_dicts_cusp,
            "ctc_param": "dt_eff",
            "idx_list": idx_list,
            "peaks_kev": peaks_kev,
        },
        {
            "peak_dicts": kwarg_dicts_zac,
            "ctc_param": "dt_eff",
            "idx_list": idx_list,
            "peaks_kev": peaks_kev,
        },
        {
            "peak_dicts": kwarg_dicts_trap,
            "ctc_param": "dt_eff",
            "idx_list": idx_list,
            "peaks_kev": peaks_kev,
        },
    ]

    fom = eval(opt_dict["fom"])
    out_field = opt_dict["fom_field"]
    out_err_field = opt_dict["fom_err_field"]
    sample_x = np.array(opt_dict["initial_samples"])

    results_cusp = []
    results_zac = []
    results_trap = []

    sample_y_cusp = []
    sample_y_zac = []
    sample_y_trap = []

    err_y_cusp = []
    err_y_zac = []
    err_y_trap = []

    for i, x in enumerate(sample_x):
        db_dict["cusp"]["sigma"] = f"{x[0]}*us"
        db_dict["zac"]["sigma"] = f"{x[0]}*us"
        db_dict["etrap"]["rise"] = f"{x[0]}*us"

        log.info(f"Initialising values {i+1} : {db_dict}")

        tb_out = run_one_dsp(tb_data, dsp_config, db_dict=db_dict, verbosity=0)

        res = fom(tb_out, kwarg_dict[0])
        results_cusp.append(res)
        sample_y_cusp.append(res[out_field])
        err_y_cusp.append(res[out_err_field])

        res = fom(tb_out, kwarg_dict[1])
        results_zac.append(res)
        sample_y_zac.append(res[out_field])
        err_y_zac.append(res[out_err_field])

        res = fom(tb_out, kwarg_dict[2])
        results_trap.append(res)
        sample_y_trap.append(res[out_field])
        err_y_trap.append(res[out_err_field])

        log.info(f"{i+1} Finished")

    if np.isnan(sample_y_cusp).all():
        max_cusp = opt_dict["nan_default"]
    else:
        max_cusp = np.ceil(np.nanmax(sample_y_cusp) * 2)
    if np.isnan(sample_y_zac).all():
        max_zac = opt_dict["nan_default"]
    else:
        max_zac = np.ceil(np.nanmax(sample_y_zac) * 2)
    if np.isnan(sample_y_trap).all():
        max_trap = opt_dict["nan_default"]
    else:
        max_trap = np.ceil(np.nanmax(sample_y_trap) * 2)

    nan_vals = [max_cusp, max_zac, max_trap]

    for i in range(len(sample_x)):
        if np.isnan(sample_y_cusp[i]):
            results_cusp[i]["y_val"] = max_cusp
            sample_y_cusp[i] = max_cusp

        if np.isnan(sample_y_zac[i]):
            results_zac[i]["y_val"] = max_zac
            sample_y_zac[i] = max_zac

        if np.isnan(sample_y_trap[i]):
            results_trap[i]["y_val"] = max_trap
            sample_y_trap[i] = max_trap

    kernel = (
        ker.ConstantKernel(2.0, constant_value_bounds="fixed")
        + 1.0 * ker.RBF(1.0, length_scale_bounds=[0.5, 2.5])
        + ker.WhiteKernel(noise_level=0.1, noise_level_bounds=(1e-5, 1e1))
    )

    lambda_param = 5
    sampling_rate = tb_data["waveform_presummed"]["dt"][0]
    sampling_unit = ureg.Quantity(tb_data["waveform_presummed"]["dt"].attrs["units"])
    waveform_sampling = sampling_rate * sampling_unit

    bopt_cusp = BayesianOptimizer(
        acq_func=opt_dict["acq_func"],
        batch_size=opt_dict["batch_size"],
        kernel=kernel,
        sampling_rate=waveform_sampling,
        fom_value=out_field,
        fom_error=out_err_field,
    )
    bopt_cusp.lambda_param = lambda_param
    bopt_cusp.add_dimension("cusp", "sigma", 0.5, 16, True, "us")

    bopt_zac = BayesianOptimizer(
        acq_func=opt_dict["acq_func"],
        batch_size=opt_dict["batch_size"],
        kernel=kernel,
        sampling_rate=waveform_sampling,
        fom_value=out_field,
        fom_error=out_err_field,
    )
    bopt_zac.lambda_param = lambda_param
    bopt_zac.add_dimension("zac", "sigma", 0.5, 16, True, "us")

    bopt_trap = BayesianOptimizer(
        acq_func=opt_dict["acq_func"],
        batch_size=opt_dict["batch_size"],
        kernel=kernel,
        sampling_rate=waveform_sampling,
        fom_value=out_field,
        fom_error=out_err_field,
    )
    bopt_trap.lambda_param = lambda_param
    bopt_trap.add_dimension("etrap", "rise", 1, 12, True, "us")

    bopt_cusp.add_initial_values(x_init=sample_x, y_init=sample_y_cusp, yerr_init=err_y_cusp)
    bopt_zac.add_initial_values(x_init=sample_x, y_init=sample_y_zac, yerr_init=err_y_zac)
    bopt_trap.add_initial_values(x_init=sample_x, y_init=sample_y_trap, yerr_init=err_y_trap)

    best_idx = np.nanargmin(sample_y_cusp)
    bopt_cusp.optimal_results = results_cusp[best_idx]
    bopt_cusp.optimal_x = sample_x[best_idx]

    best_idx = np.nanargmin(sample_y_zac)
    bopt_zac.optimal_results = results_zac[best_idx]
    bopt_zac.optimal_x = sample_x[best_idx]

    best_idx = np.nanargmin(sample_y_trap)
    bopt_trap.optimal_results = results_trap[best_idx]
    bopt_trap.optimal_x = sample_x[best_idx]

    optimisers = [bopt_cusp, bopt_zac, bopt_trap]

    out_param_dict, out_results_list = run_bayesian_optimisation(
        tb_data,
        dsp_config,
        [fom],
        optimisers,
        fom_kwargs=kwarg_dict,
        db_dict=db_dict,
        nan_val=nan_vals,
        n_iter=opt_dict["n_iter"],
    )

    Props.add_to(db_dict, out_param_dict)

    # db_dict.update(out_param_dict)

    t2 = time.time()
    log.info(f"Optimiser finished in {(t2-t1)/60} minutes")

    out_alpha_dict = {}
    out_alpha_dict["cuspEmax_ctc"] = {
        "expression": "cuspEmax*(1+dt_eff*a)",
        "parameters": {"a": float(round(bopt_cusp.optimal_results["alpha"], 9))},
    }

    out_alpha_dict["cuspEftp_ctc"] = {
        "expression": "cuspEftp*(1+dt_eff*a)",
        "parameters": {"a": float(round(bopt_cusp.optimal_results["alpha"], 9))},
    }

    out_alpha_dict["zacEmax_ctc"] = {
        "expression": "zacEmax*(1+dt_eff*a)",
        "parameters": {"a": float(round(bopt_zac.optimal_results["alpha"], 9))},
    }

    out_alpha_dict["zacEftp_ctc"] = {
        "expression": "zacEftp*(1+dt_eff*a)",
        "parameters": {"a": float(round(bopt_zac.optimal_results["alpha"], 9))},
    }

    out_alpha_dict["trapEmax_ctc"] = {
        "expression": "trapEmax*(1+dt_eff*a)",
        "parameters": {"a": float(round(bopt_trap.optimal_results["alpha"], 9))},
    }

    out_alpha_dict["trapEftp_ctc"] = {
        "expression": "trapEftp*(1+dt_eff*a)",
        "parameters": {"a": float(round(bopt_trap.optimal_results["alpha"], 9))},
    }
    if "ctc_params" in db_dict:
        db_dict["ctc_params"].update(out_alpha_dict)
    else:
        db_dict.update({"ctc_params": out_alpha_dict})

    Path(args.qbb_grid_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.qbb_grid_path).open("wb") as f:
        pkl.dump(optimisers, f)

else:
    Path(args.qbb_grid_path).touch()

Path(args.final_dsp_pars).parent.mkdir(parents=True, exist_ok=True)
Props.write_to(args.final_dsp_pars, db_dict)

if args.plot_path:
    if args.inplots:
        with Path(args.inplots).open("rb") as r:
            plot_dict = pkl.load(r)
    else:
        plot_dict = {}

    plot_dict["trap_optimisation"] = {
        "kernel_space": bopt_trap.plot(init_samples=sample_x),
        "acq_space": bopt_trap.plot_acq(init_samples=sample_x),
    }

    plot_dict["cusp_optimisation"] = {
        "kernel_space": bopt_cusp.plot(init_samples=sample_x),
        "acq_space": bopt_cusp.plot_acq(init_samples=sample_x),
    }

    plot_dict["zac_optimisation"] = {
        "kernel_space": bopt_zac.plot(init_samples=sample_x),
        "acq_space": bopt_zac.plot_acq(init_samples=sample_x),
    }

    Path(args.plot_path).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.plot_path).open("wb") as w:
        pkl.dump(plot_dict, w, protocol=pkl.HIGHEST_PROTOCOL)
