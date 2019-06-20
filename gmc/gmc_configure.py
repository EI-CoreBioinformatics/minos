import csv
import sys
import yaml
import argparse
import os
import pathlib
import subprocess

from collections import OrderedDict

from . import __version__

STRANDINFO = {"_xx", "_rf", "_fr"}

class ScoringMetricsManager(object):
	def __importMetricsData(self, fn, use_tpm=False):
		self.metrics = OrderedDict()
		for row in csv.reader(open(fn), delimiter="\t", quotechar='"'):
			if row[0].startswith("#"):
				continue
			# metric_name_prefix    metric_class    multiplier  not_fragmentary_min_value   file_path
			data = {
				"multiplier": row[2],
				"min_value": row[3],
				"files": row[4].split(",")
			}

			if row[1] == "expression":
				if use_tpm:
					if row[0][-3:] not in STRANDINFO:
						raise ValueError("ERROR: expression metric does not have strandedness information. Please add a suffix to indicate strandedness (_xx: unstranded, _rf: rf-stranded, _fr: fr-stranded). " + row[0])
				else:
					print("Found expression metric: {} but --use-tpm-for-picking was not set. Ignoring.".format(row[0]))
					continue

			self.metrics.setdefault(row[1], OrderedDict()).setdefault(row[0], list()).append(data)
	
		if len(self.metrics.get("junction", OrderedDict())) > 1:
			raise ValueError("ERROR: More than one junction file detected. " + self.metrics.get("junction", list()))

		for k in self.metrics:
			print(k)
			for j in self.metrics[k]:
				print(j, self.metrics[k][j], sep="\n")

	def __init__(self, metrics_file, scoring_template_file, outdir, prefix, use_tpm=False):
		self.__importMetricsData(metrics_file, use_tpm=use_tpm)

	def getMetricsData(self, metric):
		mdata = dict()
		for k in self.metrics.get(metric, OrderedDict()):
			mdata.setdefault(k, list()).extend(item["files"] for item in self.metrics[metric][k])
			pass
		return mdata
	
	def generateScoringFile(self, scoring_template, outfile):
		ext_metrics = list()
		for mclass in self.metrics:
			ext_metrics.extend("external.{}_aF1".format(run) for run in self.metrics[mclass])

		def generate_nf_expression(metrics):
			ext_metrics = list()
			blast_metrics = {"seq_prot", "blastdb_prot"}
			for mclass in self.metrics:
				if mclass not in {"junction", "expression"}:
					suffixes = ["qCov", "tCov"] if mclass in blast_metrics else ["aF1"]
					for run in metrics[mclass]:
						for suffix in suffixes:						
							ext_metrics.append("external.{}_{}".format(run, suffix))				
					#ext_metrics.extend("external.{}_aF1".format(run) for run in metrics[mclass])

			expression = "[((exon_num.multi and (combined_cds_length.multi or {0}))"
			expression += ", or, "
			expression += "(exon_num.mono and (combined_cds_length.mono or {0})))]"
			return expression.format("*".join(ext_metrics).replace("*", " or "))

		def generate_nf_params(metrics):
			params = list() 
			blast_metrics = {"seq_prot", "blastdb_prot"}
			for mclass in metrics:
				if mclass not in {"junction", "expression"}:
					suffixes = ["qCov", "tCov"] if mclass in blast_metrics else ["aF1"]
					for run in metrics[mclass]:
						for suffix in suffixes:
							params.append("    external.{}_{}: {{operator: gt, value: {}}}".format(run, suffix, metrics[mclass][run][0]["min_value"]))

			return params
			
		def generate_external_scoring(metrics):
			blast_metrics = {"seq_prot", "blastdb_prot"}
			scoring = list()			

			for mclass in metrics:
				if mclass == "expression":
					suffixes = ["tpm"]
				elif mclass in blast_metrics:
					suffixes = ["qCov", "tCov"]
				elif mclass == "junction":
					continue
				else:
					suffixes = ["nF1", "jF1", "eF1", "aF1"]					
					
				for run in metrics[mclass]:
					for suffix in suffixes:
						comment = "#" if suffix in ["nF1", "jF1", "eF1"] else ""
						if suffix == "tpm":
							expression = "{{rescaling: min, filter: {{operator: lt, value: {}}}, multiplier: {}}}".format(
								0.5, #! GET VALUE!
								metrics[mclass][run][0]["multiplier"]
							)
							suffix = ""
						else:
							expression = "{{rescaling: max, use_raw: true, multiplier: {}}}".format(metrics[mclass][run][0]["multiplier"])
							suffix = "_" + suffix
						scoring.append(
							"  {}external.{}{}: {}".format(comment, run, suffix, expression)
						)
			
			return scoring		
				
		

		with open(scoring_template) as _in, open(outfile, "wt") as _out:
			for line in _in:
				if line.strip().startswith("### GMC_CONFIGURE ###"):
					break
				print(line, end="", file=_out)

			print("not_fragmentary:", file=_out)
			print("  # expression: [combined_cds_length]", file=_out)
			print("  expression:", generate_nf_expression(self.metrics), file=_out)
			print("  parameters:", file=_out)
			print("    # is_complete: {operator: eq, value: true}", file=_out)
			print("    exon_num.multi: {operator: gt, value: 1}", file=_out)
			print("    # cdna_length.multi: {operator: ge, value: 200}", file=_out)
			print("    combined_cds_length.multi: {operator: gt, value: 200}", file=_out)
			print("    exon_num.mono: {operator: eq, value: 1}", file=_out)              			
			print("    combined_cds_length.mono: {operator: gt, value: 300}", file=_out) 
			print("    # combined_cds_length: {operator: gt, value: 300}", file=_out)
			# print("    external.all_aF1: {operator: gt, value: 0.5}", file=_out)
			print(*generate_nf_params(self.metrics), sep="\n", file=_out)
			print("scoring:", file=_out)
			print("  # external metrics START", file=_out)
			print(*generate_external_scoring(self.metrics), sep="\n", file=_out)
			# print("  # external.tpsi_cov: {rescaling: max, use_raw: true, multiplier: 10}", file=_out)
			# print("  # external.all_repeats_cov: {rescaling: max, use_raw: true, multiplier: 10}", file=_out)
			# print("  # external.interspersed_repeats_cov: {rescaling: max, use_raw: true, multiplier: 10}", file=_out)
			"""
			external.EI_tpm: {rescaling: min, filter: {operator: lt, value: 0.5}, multiplier: 1}
			So as discussed this replaces the following lines
			  # external.EI_tpm_05: {rescaling: max, use_raw: true, multiplier: 10}
			  # external.EI_tpm_1: {rescaling: max, use_raw: true, multiplier: 10}
			"""
			print("  external.cpc: {rescaling: max, use_raw: true, multiplier: 1}", file=_out)
			# print("  # all boolean metrics values from here below", file=_out)
			# print("  # external.EI_tpm_05: {rescaling: max, use_raw: true, multiplier: 10}", file=_out)
			# print("  # external.EI_tpm_1: {rescaling: max, use_raw: true, multiplier: 10}", file=_out)
			# print("  # external.SRA_tpm_05: {rescaling: max, use_raw: true, multiplier: 10}", file=_out)
			# print("  # external.SRA_tpm_1: {rescaling: max, use_raw: true, multiplier: 10}", file=_out)
			# print("  external.fln: {rescaling: max, use_raw: true, multiplier: 5}", file=_out)
			print("  # external metrics END", file=_out)
			
			for line in _in:
				print(line, end="", file=_out)

	pass



def run_configure(args):

	pathlib.Path(args.outdir).mkdir(exist_ok=True, parents=True)
	pathlib.Path(os.path.join(args.outdir, "hpc_logs")).mkdir(exist_ok=True, parents=True)

	scoringFile = os.path.join(args.outdir, args.prefix + ".scoring.yaml")
	smm = ScoringMetricsManager(args.external_metrics, args.scoring_template, args.outdir, args.prefix, use_tpm=args.use_tpm_for_picking)
	print("Generating scoring file " + scoringFile + " ...", end="", flush=True)
	smm.generateScoringFile(args.scoring_template, scoringFile)
	print(" done.")

	#!TODO: 
	# - scan config template for reference
	# - warn if not present
	# - add command line option 
	
	mikado_config_file = os.path.join(args.outdir, args.prefix + ".mikado_config.yaml")


	cmd = "singularity exec {} mikado configure --list {}{}-od {} --scoring {}{}{}".format(
		args.mikado_container,
		args.list_file,
		(" --external " + args.external + " ") if args.external else " ",
		args.outdir,
		scoringFile,
		(" --junctions " + list(smm.getMetricsData("junction").values())[0][0][0] + " ") if smm.getMetricsData("junction") else " ",
		mikado_config_file
	)

	print(cmd)
	out = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)

	run_zzz_config = {
		"prefix": args.prefix,
		"outdir": args.outdir,
		"mikado-container": args.mikado_container,
		"mikado-config-file": mikado_config_file,
		"external-metrics": args.external_metrics,
		"reference-sequence": args.reference,
		"blast-mode": args.blastmode,
		"use-tpm-for-picking": args.use_tpm_for_picking,
		"external-metrics-data": args.external_metrics,
		"annotation_version": args.annotation_version,
		"genus_identifier": args.genus_identifier
	}

	run_zzz_data = { 
		"data": {
			"expression-runs": smm.getMetricsData("expression"), 
			"transcript-runs": smm.getMetricsData("aln_tran"), 
			"protein-runs": smm.getMetricsData("aln_prot"), 
			"protein-seqs": smm.getMetricsData("seq_prot"),
			"junction-data": smm.getMetricsData("junction"),	
			"repeat-data": smm.getMetricsData("repeat")
		}
	}


	
	with open(os.path.join(args.outdir, args.prefix + ".run_config.yaml"), "wt") as run_config_out:
		yaml.dump(run_zzz_config, run_config_out, default_flow_style=False)
		yaml.dump(run_zzz_data, run_config_out, default_flow_style=False)

	pass

