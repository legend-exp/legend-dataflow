import argparse
import pickle as pkl
import time
from pathlib import Path

import lgdo.lh5 as lh5
import numpy as np
import pygama.pargen.noise_optimization as pno
from dbetto import TextDB
from dbetto.catalog import Props
from pygama.pargen.data_cleaning import generate_cuts, get_cut_indexes
from pygama.pargen.dsp_optimize import run_one_dsp

from .....log import build_log


def par_geds_dsp_nopt() -> None:
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--raw-filelist", help="raw_filelist", type=str)
    argparser.add_argument("--database", help="database", type=str, required=True)
    argparser.add_argument("--inplots", help="inplots", type=str)

    argparser.add_argument("--configs", help="configs", type=str, required=True)
    argparser.add_argument("--log", help="log_file", type=str)

    argparser.add_argument("--datatype", help="Datatype", type=str, required=True)
    argparser.add_argument("--timestamp", help="Timestamp", type=str, required=True)
    argparser.add_argument("--channel", help="Channel", type=str, required=True)
    argparser.add_argument(
        "--raw-table-name", help="raw table name", type=str, required=True
    )

    argparser.add_argument("--dsp-pars", help="dsp_pars", type=str, required=True)
    argparser.add_argument("--plot-path", help="plot_path", type=str)

    args = argparser.parse_args()

    configs = TextDB(args.configs, lazy=True).on(args.timestamp, system=args.datatype)
    config_dict = configs["snakemake_rules"]["pars_dsp_nopt"]

    log = build_log(config_dict, args.log)

    t0 = time.time()

    dsp_config = config_dict["inputs"]["processing_chain"][args.channel]
    opt_json = config_dict["inputs"]["optimiser_config"][args.channel]

    opt_dict = Props.read_from(opt_json)
    db_dict = Props.read_from(args.database)

    if opt_dict.pop("run_nopt") is True:
        with Path(args.raw_filelist).open() as f:
            files = f.read().splitlines()

        raw_files = sorted(files)

        energies = lh5.read_as(
            f"{args.raw_table_name}/daqenergy", raw_files, library="np"
        )
        idxs = np.where(energies == 0)[0]
        tb_data = lh5.read(
            args.raw_table_name, raw_files, n_rows=opt_dict["n_events"], idx=idxs
        )
        t1 = time.time()
        log.info(f"Time to open raw files {t1 - t0:.2f} s, n. baselines {len(tb_data)}")

        log.info(f"Select baselines {len(tb_data)}")
        dsp_data = run_one_dsp(tb_data, dsp_config)
        cut_dict = generate_cuts(dsp_data, cut_dict=opt_dict.pop("cut_pars"))
        cut_idxs = get_cut_indexes(dsp_data, cut_dict)
        tb_data = lh5.read(
            args.raw_table_name,
            raw_files,
            n_rows=opt_dict.pop("n_events"),
            idx=idxs[cut_idxs],
        )
        log.info(f"... {len(tb_data)} baselines after cuts")

        if isinstance(dsp_config, (str, list)):
            dsp_config = Props.read_from(dsp_config)

        if args.plot_path:
            out_dict, plot_dict = pno.noise_optimization(
                tb_data,
                dsp_config,
                db_dict.copy(),
                opt_dict,
                args.raw_table_name,
                display=1,
            )
        else:
            out_dict = pno.noise_optimization(
                raw_files, dsp_config, db_dict.copy(), opt_dict, args.raw_table_name
            )

        t2 = time.time()
        log.info(f"Optimiser finished in {(t2 - t0) / 60} minutes")
    else:
        out_dict = {}
        plot_dict = {}

    if args.plot_path:
        Path(args.plot_path).parent.mkdir(parents=True, exist_ok=True)
        if args.inplots:
            with Path(args.inplots).open("rb") as r:
                old_plot_dict = pkl.load(r)
            plot_dict = dict(noise_optimisation=plot_dict, **old_plot_dict)
        else:
            plot_dict = {"noise_optimisation": plot_dict}
        with Path(args.plot_path).open("wb") as f:
            pkl.dump(plot_dict, f, protocol=pkl.HIGHEST_PROTOCOL)

    Path(args.dsp_pars).parent.mkdir(parents=True, exist_ok=True)
    Props.write_to(args.dsp_pars, dict(nopt_pars=out_dict, **db_dict))
