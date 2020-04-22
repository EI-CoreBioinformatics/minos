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

from gmc import __version__
from gmc.gmc_configure import *
from eicore.snakemake_helper import *


def add_default_options(parser):
	common_group = parser.add_argument_group("gmc options")
			
	# common_group.add_argument("list_file", type=str)
	common_group.add_argument("--outdir", "-o", type=str, default="gmc_run")
	common_group.add_argument("--prefix", type=str, default="gmc_run")
	common_group.add_argument("--mikado-container", type=str, default="/ei/software/testing/gmc/dev/x86_64/mikado.simg")
	#common_group.add_argument("--mikado-container", type=str, default="/ei/software/testing/mikado/20190325_c940de1/x86_64/mikado-20190325_c940de1.simg")
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
	configure_parser.add_argument("--annotation-version", type=str, default="EIv1")
	configure_parser.add_argument("--genus-identifier", type=str, default="XYZ")
	configure_parser.add_argument("--use-tpm-for-picking", action="store_true")
	configure_parser.add_argument("--force-reconfiguration", "-f", action="store_true")
	configure_parser.add_argument("--config-file", type=str, default=DEFAULT_CONFIG_FILE)
	configure_parser.add_argument("--busco-level", type=str, default="all") # proteins in release 
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

	run_configuration_file = None
	try:
		run_configuration_file = os.path.abspath(glob.glob(os.path.join(args.outdir, "*.run_config.yaml")).pop())
	except:
		pass 


	print(args.runmode)
	if args.runmode == "configure":
		if run_configuration_file is None or args.force_reconfiguration:
			run_configure(args)
		elif run_configuration_file is not None:
			print("Configuration file {} already present. Please set --force-reconfiguration/-f to override this.".format(run_configuration_file))
	elif args.runmode == "run":
		snake = join(dirname(__file__), "zzz", "gmc_run.smk")
		if run_configuration_file is None:
			raise ValueError("Missing run configuration in " + args.outdir)

		print("Using run configuration file: " + run_configuration_file)

		if args.mikado_container and os.path.exists(args.mikado_container):
			rconf = yaml.load(open(run_configuration_file), Loader=yaml.SafeLoader)
			rconf["mikado-container"] = args.mikado_container
			
			run_configuration_file = run_configuration_file.replace(".yaml", ".with_singularity.yaml")
			with open(run_configuration_file, "wt") as run_config_out:
				yaml.dump(rconf, run_config_out, default_flow_style=False)


		exe_env = ExecutionEnvironment(args, NOW, job_suffix="GMC_" + args.outdir, log_dir=os.path.join(args.outdir, "hpc_logs"))

		run_snakemake(snake, args.outdir, run_configuration_file, exe_env, dryrun=args.dryrun)
	

	pass

if __name__ == "__main__": 
	main()
