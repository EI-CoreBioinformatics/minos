import csv
import sys
import yaml
import argparse
import os
from os.path import join, basename, dirname, abspath
import pathlib
import subprocess
import glob
import shutil

from collections import OrderedDict

from minos import __version__
from minos.minos_configure import *
from eicore.snakemake_helper import *


def add_default_options(parser):
	common_group = parser.add_argument_group("minos options")
			
	common_group.add_argument("--outdir", "-o", type=str, default="minos_run")
	common_group.add_argument("--prefix", type=str, default="minos_run")
	common_group.add_argument("--mikado-container", type=str, default="/ei/software/testing/minos/dev/x86_64/mikado.simg")
	common_group.add_argument("--dryrun", action="store_true")
	make_exeenv_arg_group(parser, allow_mode_selection=False, silent=True)
		
def add_configure_parser(subparsers):
	configure_parser = subparsers.add_parser(
		"configure",
		help="",
		description=""
	)
	
	configure_parser.add_argument("list_file", type=str)	
	configure_parser.add_argument("scoring_template", type=str)
	configure_parser.add_argument("reference", type=str)
	configure_parser.add_argument("--external", type=str, default="")
	configure_parser.add_argument("--external-metrics", type=str, default="")
	configure_parser.add_argument("--blastmode", choices=("blastp", "blastx"), default="blastp")
	configure_parser.add_argument("--annotation-version", type=str, default="EIv1")
	configure_parser.add_argument("--genus-identifier", type=str, default="XYZ")
	configure_parser.add_argument("--use-tpm-for-picking", action="store_true")
	configure_parser.add_argument("--force-reconfiguration", "-f", action="store_true")
	configure_parser.add_argument("--config-file", type=str, default=DEFAULT_CONFIG_FILE)
	configure_parser.add_argument("--busco-level", type=str, default="all") # proteins in release 
	configure_parser.add_argument("--busco-scoring", type=int, help="Force busco protein runs and use results in transcript scoring with the specified multiplier.")
	configure_parser.add_argument("--busco-lineage", type=str, help="Required if --busco-level is not in {none,off}.")
	configure_parser.add_argument("--busco-genome-run",type=str, help="Directory with short_summary.txt and full_table.tsv from processing the reference with busco genome.")
	
	add_default_options(configure_parser)
	configure_parser.set_defaults(runmode="configure")


def add_run_parser(subparsers):
	run_parser = subparsers.add_parser(
		"run",
		help="",
		description=""
	)

	run_parser.add_argument("--rerun-from", type=str, choices=("start", "pick", "collapse_metrics", "off"), default="off", help="Rerun from specific stage [off]")

	add_default_options(run_parser)
	run_parser.set_defaults(runmode="run")

def parse_args():
	ap = argparse.ArgumentParser(prog="minos", description="The Earlham Institute Gene Model Consolidation Pipeline (minos).")
	subparsers = ap.add_subparsers(
		help=""
	)

	add_configure_parser(subparsers)
	add_run_parser(subparsers)

	return ap.parse_args()


def main():	
	print("Starting MINOS V " + __version__)
	print()
	
	if len(sys.argv) == 1:
		sys.argv.append("-h")
	
	args = parse_args()

	run_configuration_file = None
	try:
		run_configuration_file = os.path.abspath(glob.glob(os.path.join(args.outdir, "*.run_config.yaml")).pop())
	except:
		pass 

	print("Runmode is", args.runmode)
	if args.runmode == "configure":
		if run_configuration_file is None or args.force_reconfiguration:
			MinosRunConfiguration(args).run()
		elif run_configuration_file is not None:
			print("Configuration file {} already present. Please set --force-reconfiguration/-f to override this.".format(run_configuration_file))
	elif args.runmode == "run":
		snake = join(dirname(__file__), "zzz", "minos_run.smk")

		if run_configuration_file is None:
			raise ValueError("Missing run configuration in " + args.outdir)

		print("Using run configuration file: " + run_configuration_file)

		run_config = yaml.load(open(run_configuration_file), Loader=yaml.SafeLoader)
		if args.mikado_container and os.path.exists(args.mikado_container):
			run_config["mikado-container"] = args.mikado_container
		if args.hpc_config and os.path.exists(args.hpc_config) and args.hpc_config != DEFAULT_HPC_CONFIG_FILE:
			run_config["hpc_config"] = args.hpc_config
		args.hpc_config = run_config["hpc_config"]

		run_configuration_file = run_configuration_file.replace(".yaml", ".{}.yaml".format(NOW))
		with open(run_configuration_file, "wt") as run_config_out:
			yaml.dump(run_config, run_config_out, default_flow_style=False, sort_keys=False)

		exe_env = ExecutionEnvironment(args, NOW, job_suffix="MINOS_" + args.outdir, log_dir=os.path.join(args.outdir, "hpc_logs"))

		minos_complete_sentinel = os.path.join(args.outdir, "MINOS_RUN_COMPLETE")
		results_dir = os.path.join(args.outdir, "results")
		if os.path.exists(minos_complete_sentinel) and args.rerun_from != "off":

			if os.path.exists(results_dir):

				archive_dir = os.path.join(args.outdir, "archives")
				pathlib.Path(archive_dir).mkdir(exist_ok=True, parents=True)
				dest_dir = os.path.join(archive_dir, "{}_{}".format(os.path.basename(args.outdir), NOW))
				try:
					shutil.copytree(results_dir, dest_dir)
				except:
					raise ValueError("Rerunning targeting a previously completed run. {results_dir} present but could not archive it in {archive_dir}. Please (re)move results dir manually before proceeding".format(results_dir=results_dir, archive_dir=dest_dir))

			serialise_sentinel = os.path.join(args.outdir, "MIKADO_SERIALISE_DONE")
			collapse_sentinel = os.path.join(args.outdir, "COLLAPSE_METRICS_DONE")
			if args.rerun_from == "start" and os.path.exists(run_config["mikado-config-file"]):
				open(run_config["mikado-config-file"], "a").close()
			elif args.rerun_from == "pick" and os.path.exists(serialise_sentinel):
				open(serialise_sentinel, "w").close()
			elif args.rerun_from == "collapse_metrics" and os.path.exists(collapse_sentinel):
				open(collapse_sentinel, "w").close()

		result = run_snakemake(snake, args.outdir, run_configuration_file, exe_env, dryrun=args.dryrun)
		if result:
			open(minos_complete_sentinel, "w").close()
	

	pass

if __name__ == "__main__": 
	main()
