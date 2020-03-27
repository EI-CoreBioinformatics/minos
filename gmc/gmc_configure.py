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
MIKADO_CONFIGURE_CMD = "singularity exec {container} mikado configure --list {list_file}{external_metrics}-od {output_dir} --reference {reference} --scoring {scoring_file}{junctions}{mikado_config_file} --full"

class ScoringMetricsManager(object):
	def __importMetricsData(self, fn, use_tpm=False):
		self.metrics = OrderedDict()
		for row in csv.reader(open(fn), delimiter="\t", quotechar='"'):
			try:
				metric_name, metric_class, multiplier, nf_minval, files = row
			except:
				raise ValueError("External metrics record has to comply to the format: metric_name_prefix    metric_class    multiplier  not_fragmentary_min_value   file_path\n{}".format(row))

			if metric_name.startswith("#"):
				continue
			if "." in metric_name:
				raise ValueError("Metric name {} contains invalid character '.'. Please adjust metric names accordingly.".format(metric_name))
			files = files.split(",")
			for f in files:
				if not os.path.exists(f):
					raise ValueError("Missing input file for metric {}: {}".format(metric_name, f))

			data = {
				"multiplier": multiplier,
				"min_value": nf_minval,
				"files": files
			}

			if metric_class == "expression":
				if use_tpm:
					if metric_name[-3:] not in STRANDINFO:
						raise ValueError("Expression metric does not have strandedness information. Please add a suffix to indicate strandedness (_xx: unstranded, _rf: rf-stranded, _fr: fr-stranded). " + metric_name)
				else:
					print("Found expression metric: {} but --use-tpm-for-picking was not set. Ignoring.".format(metric_name))
					continue

			self.metrics.setdefault(metric_class, OrderedDict()).setdefault(metric_name, list()).append(data)
	
		if len(self.metrics.get("junction", OrderedDict())) > 1:
			raise ValueError("More than one junction file detected. " + self.metrics.get("junction", list()))

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
				if line.strip().startswith("### GMC:GENERATE_NF_EXPRESSION"):
					line = next(_in)
					print("  expression:", generate_nf_expression(self.metrics), file=_out)
				elif line.strip().startswith("### GMC:GENERATE_NF_PARAMS"):
					line = next(_in)
					# print("    external.all_aF1: {operator: gt, value: 0.5}", file=_out)
					print(*generate_nf_params(self.metrics), sep="\n", file=_out)
				elif line.strip().startswith("### GMC:GENERATE_EXTERNAL_SCORING"):
					line = next(_in)
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
				else:
					print(line, end="", file=_out)				
			
	pass



def run_configure(args):

	pathlib.Path(args.outdir).mkdir(exist_ok=True, parents=True)
	pathlib.Path(os.path.join(args.outdir, "hpc_logs")).mkdir(exist_ok=True, parents=True)

	scoring_file = os.path.join(args.outdir, args.prefix + ".scoring.yaml")
	smm = ScoringMetricsManager(args.external_metrics, args.scoring_template, args.outdir, args.prefix, use_tpm=args.use_tpm_for_picking)
	print("Generating scoring file " + scoring_file + " ...", end="", flush=True)
	smm.generateScoringFile(args.scoring_template, scoring_file)
	print(" done.")

	#!TODO: 
	# - scan config template for reference
	# - warn if not present
	# - add command line option 
	
	mikado_config_file = os.path.join(args.outdir, args.prefix + ".mikado_config.yaml")
	gmc_config = yaml.load(open(args.config_file), Loader=yaml.SafeLoader)

	cmd = MIKADO_CONFIGURE_CMD.format(
		container=args.mikado_container,
		list_file=args.list_file,
		external_metrics=(" --external " + args.external + " ") if args.external else " ",
		output_dir=args.outdir,
		reference=args.reference,
		scoring_file=scoring_file,
		junctions=(" --junctions " + list(smm.getMetricsData("junction").values())[0][0][0] + " ") if smm.getMetricsData("junction") else " ",
		mikado_config_file=mikado_config_file
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
		yaml.dump(gmc_config, run_config_out, default_flow_style=False)

	pass

