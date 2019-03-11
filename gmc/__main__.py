import csv
import sys
import yaml
import argparse
import os
from os.path import join, basename, dirname, abspath
import pathlib
import subprocess
import glob

from collections import OrderedDict

from . import __version__
from .gmc_configure import *
from .snakemake_helper import *


def add_default_options(parser):
	common_group = parser.add_argument_group("gmc options")
			
	# common_group.add_argument("list_file", type=str)
	common_group.add_argument("--outdir", "-o", type=str, default="gmc_run")
	common_group.add_argument("--prefix", type=str, default="gmc_run")
	common_group.add_argument("--mikado-container", type=str, default="/ei/software/testing/gmc/dev/x86_64/mikado.simg")
	common_group.add_argument("--dryrun", action="store_true")
	make_exeenv_arg_group(parser, default_hpc_config_file="", allow_mode_selection=False, silent=True)
		
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
	
	add_default_options(configure_parser)
	configure_parser.set_defaults(runmode="configure")


def add_run_parser(subparsers):
	run_parser = subparsers.add_parser(
		"run",
		help="",
		description=""
	)

	add_default_options(run_parser)
	run_parser.set_defaults(runmode="run")
		



def main():	
	print("Starting EI GMC V " + __version__)
	print()
	
	if len(sys.argv) == 1:
		sys.argv.append("-h")
	
	ap = argparse.ArgumentParser(prog="gmc", description="The Earlham Institute Gene Model Consolidation Pipeline (gmc).")
	
	subparsers = ap.add_subparsers(
		help=""
	)

	add_configure_parser(subparsers)
	add_run_parser(subparsers)

	args = ap.parse_args()


	print(args.runmode)
	if args.runmode == "configure":
		run_configure(args)	
	elif args.runmode == "run":
		snake = join(dirname(__file__), "zzz", "gmc_run.smk.py")
		try:
			run_config = abspath(glob.glob(os.path.join(args.outdir, "*.run_config.yaml")).pop())
		except:
			raise ValueError("Missing run configuration in " + args.outdir)

		print("Using run configuration file: " + run_config)

		exe_env = ExecutionEnvironment(args, NOW, job_suffix="GMC_" + args.outdir, log_dir=os.path.join(args.outdir, "hpc_logs"))

		run_snakemake(snake, args.outdir, run_config, exe_env, dryrun=args.dryrun)
	

	pass

if __name__ == "__main__": 
	main()
