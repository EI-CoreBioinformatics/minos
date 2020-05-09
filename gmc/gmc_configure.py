import csv
import sys
import yaml
import os
import glob
import pathlib
import subprocess

from enum import Enum, unique, auto

from gmc import __version__
from gmc.gmc_scoring import ScoringMetricsManager
from gmc.busco_configure import BuscoConfiguration

#!TODO: 
# - scan config template for reference
# - warn if not present
# - add command line option 


@unique
class ExternalMetrics(Enum):
    MIKADO_TRANSCRIPTS_OR_PROTEINS = auto()
    PROTEIN_BLAST_TOPHITS = auto()
    CPC_CODING_POTENTIAL = auto()
    KALLISTO_TPM_EXPRESSION = auto()
    REPEAT_ANNOTATION = auto()
    BUSCO_PROTEINS = auto()

MIKADO_CONFIGURE_CMD = "{cmd} --list {list_file}{external_metrics}-od {output_dir} --reference {reference} --scoring {scoring_file}{junctions}{mikado_config_file} --full"


class GmcRunConfiguration(dict):
	def _run_mikado_configure(self, args):
		cmd = MIKADO_CONFIGURE_CMD.format(
			cmd=self["program_calls"]["mikado"].format(container=args.mikado_container, program="configure"),
			list_file=args.list_file,
			external_metrics=(" --external " + args.external + " ") if args.external else " ",
			output_dir=args.outdir,
			reference=args.reference,
			scoring_file=self.scoring_file,
			junctions=(" --junctions " + list(self.smm.getMetricsData("junction").values())[0][0][0] + " ") if self.smm.getMetricsData("junction") else " ",
			mikado_config_file=self.mikado_config_file
		)
		print("Running mikado configure with:", cmd, sep="\n")
		out = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
		print(out)
	
	def _generate_scoring_file(self, args):
		self.scoring_file = os.path.join(args.outdir, args.prefix + ".scoring.yaml")
		self.smm = ScoringMetricsManager(args)
		print("Generating scoring file " + self.scoring_file + " ...", end="", flush=True)
		self.smm.generateScoringFile(args.scoring_template, self.scoring_file, busco_scoring=args.busco_scoring)
		print(" done.")                                                                                                                      	
	def __init__(self, args):
		print("Configuring run...")
		print(args)
		self.args = args
		pathlib.Path(args.outdir).mkdir(exist_ok=True, parents=True)
		pathlib.Path(os.path.join(args.outdir, "hpc_logs")).mkdir(exist_ok=True, parents=True)
		self.mikado_config_file = os.path.join(args.outdir, args.prefix + ".mikado_config.yaml")



	def _generate_run_configuration(self, args):
		print("Generating run configuration file ...", flush=True, end="")
		self.update({
			"prefix": args.prefix,
			"outdir": args.outdir,
			"mikado-container": args.mikado_container,
			"mikado-config-file": self.mikado_config_file,
			"external-metrics": args.external_metrics,
			"reference-sequence": args.reference,
			"blast-mode": args.blastmode,
			"use-tpm-for-picking": args.use_tpm_for_picking,
			"external-metrics-data": args.external_metrics,
			"annotation_version": args.annotation_version,
			"genus_identifier": args.genus_identifier,
			"hpc_config": args.hpc_config
		})
		self["busco_analyses"] = dict(BuscoConfiguration(args).items())

		self["data"] = {"transcript_models": {row[1]: row[0] for row in csv.reader(open(args.list_file), delimiter="\t")}}
		self["data"].update(self.smm.get_scoring_data())
		self.update(yaml.load(open(args.config_file), Loader=yaml.SafeLoader))

		with open(os.path.join(args.outdir, args.prefix + ".run_config.yaml"), "wt") as run_config_out:
			yaml.dump(dict(self.items()), run_config_out, default_flow_style=False, sort_keys=False)
		print(" done.")

	def run(self):
		self._generate_scoring_file(self.args)
		self._generate_run_configuration(self.args)
		self._run_mikado_configure(self.args)
