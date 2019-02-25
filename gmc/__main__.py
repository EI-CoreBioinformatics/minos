import csv
import sys
import yaml
import argparse
import os
import pathlib
import subprocess

from collections import OrderedDict

from . import __version__
from .gmc_configure import *
from .snakemake_helper import *


def add_default_options(parser):
	common_group = parser.add_argument_group("gmc options")
			
	common_group.add_argument("list_file", type=str)
	common_group.add_argument("--outdir", "-o", type=str, default="gmc_run")
	common_group.add_argument("--prefix", type=str, default="gmc_run")
	common_group.add_argument("--mikado-container", type=str, default="/ei/software/testing/gmc/dev/x86_64/mikado.simg")
	make_exeenv_arg_group(parser, default_hpc_config_file="", allow_mode_selection=False, silent=True)
		
def add_configure_parser(subparsers):
	configure_parser = subparsers.add_parser(
		"configure",
		help="",
		description=""
	)
	
	configure_parser.add_argument("scoring_template", type=str)
	configure_parser.add_argument("--external", type=str, default="")

	
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


	"""
	--external EXTERNAL   External configuration file to overwrite/add values from.
                            Parameters specified on the command line will take precedence over those present in the configuration file.	
	--list LIST           List of the inputs, one by line, in the form:
                        <file1>  <label>  <strandedness (true/false)> <bonus/malus (default 0)>
	 --reference REFERENCE
                        Fasta genomic reference.
	-od OUT_DIR, --out-dir OUT_DIR
                        Destination directory for the output.
	--scoring SCORING


	-i INTRON_RANGE INTRON_RANGE, --intron-range INTRON_RANGE INTRON_RANGE
                        Range into which intron lengths should fall, as a couple of integers.
                                                     Transcripts with intron lengths outside of this range will be penalised.
                                                     Default: (60, 900)
	  --pad                 Whether to pad transcripts in loci.

	--junctions JUNCTIONS
	  -bt BLAST_TARGETS, --blast_targets BLAST_TARGETS

  --daijin              Flag. If set, the configuration file will be also valid for Daijin.
  -bc BLAST_CHUNKS, --blast-chunks BLAST_CHUNKS
                        Number of parallel DIAMOND/BLAST jobs to run. Default: 10.
  --use-blast           Flag. If switched on, Mikado will use BLAST instead of DIAMOND.
  --use-transdecoder    Flag. If switched on, Mikado will use TransDecoder instead of Prodigal.
  --mode {nosplit,stringent,lenient,permissive,split} [{nosplit,stringent,lenient,permissive,split} ...]
                        Mode(s) in which Mikado will treat transcripts with multiple ORFs.
                        - nosplit: keep the transcripts whole.
                        - stringent: split multi-orf transcripts if two consecutive ORFs have both BLAST hits
                                     and none of those hits is against the same target.
                        - lenient: split multi-orf transcripts as in stringent, and additionally, also when
                                   either of the ORFs lacks a BLAST hit (but not both).
                        - permissive: like lenient, but also split when both ORFs lack BLAST hits
                        - split: split multi-orf transcripts regardless of what BLAST data is available.
                        If multiple modes are specified, Mikado will create a Daijin-compatible configuration file.
	"""
	args = ap.parse_args()


	print(args.runmode)
	if args.runmode == "configure":
		run_configure(args)	
	"""
	listFile = parseListFile(args.list_file)

	pathlib.Path(args.outdir).mkdir(exist_ok=True, parents=True)
	scoringFile = os.path.join(args.outdir, args.prefix + ".scoring.yaml")

	createScoringFile(args.scoring_template, listFile, scoringFile)

	" mikado configure --list list.txt --reference chr5.fas --mode permissive --scoring plants.yaml  --copy-scoring plants.yaml --junctions junctions.bed -bt uniprot_sprot_plants.fasta configuration.yaml"

	cmd = "singularity exec {} mikado configure --list {} {} -od {} --scoring {} {}".format(
		args.mikado_container,
		args.list_file,
		("--external " + args.external) if args.external else "",
		args.outdir,
		scoringFile,
		os.path.abspath("mikado_config.yaml")
	)

	print(cmd)
	out = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)


	with open(os.path.join(args.outdir, args.prefix + ".config.yaml"), "wt") as config_out:
		print(out.decode(), sep="\n", file=config_out)
	"""

	pass

if __name__ == "__main__": 
	main()
