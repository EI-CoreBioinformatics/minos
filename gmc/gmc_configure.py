import csv
import sys
import yaml
import argparse
import os
import glob
import pathlib
import subprocess

from collections import OrderedDict

from gmc import __version__

STRANDINFO = {"_xx": "unstranded", "_rf": "rf-stranded", "_fr": "fr-stranded"}

MIKADO_CONFIGURE_CMD = "{cmd} --list {list_file}{external_metrics}-od {output_dir} --reference {reference} --scoring {scoring_file}{junctions}{mikado_config_file} --full"

BUSCO_LEVELS = {"proteins", "proteome", "transcripts", "transcriptome", "genome", "none", "off", "p", "t", "g", "a", "all", "prot", "tran", "geno"}

EXTERNAL_METRICS_HEADERS = ["metric_name_prefix", "metric_class", "multiplier", "not_fragmentary_min_value", "file_path"]

NF_SUFFIXES = {
	"seq_prot": ["qCov", "tCov"],
	"blastdb_prot": ["qCov", "tCov"],
	"repeat": ["cov"],
	"expression": ["tpm"]
}

class ScoringMetricsManager(object):
	def __importMetricsData(self, fn, use_tpm=False):
		self.metrics = OrderedDict()
		for row in csv.reader(open(fn), delimiter="\t", quotechar='"'):
			try:
				metric_name, metric_class, multiplier, nf_minval, files = row
			except:
				raise ValueError("External metrics record has to comply to the format: {}\n{}".format(EXTERNAL_METRICS_HEADERS, row))

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
						raise ValueError("Expression metric {} does not have strandedness information. Please add a suffix to indicate strandedness ({})".format(
							metric_name, STRANDINFO
						))
				else:
					print("Found expression metric: {} but --use-tpm-for-picking was not set. Ignoring.".format(metric_name))
					continue

			self.metrics.setdefault(metric_class, OrderedDict()).setdefault(metric_name, list()).append(data)
	
		if len(self.metrics.get("junction", OrderedDict())) > 1:
			raise ValueError("More than one junction file supplied. " + self.metrics.get("junction", list()))

		#for k in self.metrics:
		#	print(k)
		#	for j in self.metrics[k]:
		#		print(j, self.metrics[k][j], sep="\n")

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
		for mclass, runs in self.metrics.items():
			ext_metrics.extend("external.{}_aF1".format(runid) for runid, run in runs.items())

		def generate_nf_expression(metrics):
			ext_metrics = list()
			for mclass, runs in self.metrics.items():
				if mclass not in {"junction", "expression", "repeat"}:
					suffixes = NF_SUFFIXES.get(mclass, ["aF1"])
					for runid, run in runs.items():
						for suffix in suffixes:						
							ext_metrics.append("external.{}_{}".format(runid, suffix))				

			expression = "[((exon_num.multi and (combined_cds_length.multi or {0}))"
			expression += ", or, "
			expression += "(exon_num.mono and (combined_cds_length.mono or {0})))]"
			return expression.format("*".join(ext_metrics).replace("*", " or "))

		def generate_nf_params(metrics):
			params = list() 
			for mclass, runs in metrics.items():
				if mclass not in {"junction", "expression", "repeat"}:
					suffixes = NF_SUFFIXES.get(mclass, ["aF1"])
					for runid, run in runs.items():
						for suffix in suffixes:
							params.append("    external.{}_{}: {{operator: gt, value: {}}}".format(runid, suffix, run[0]["min_value"]))

			return params
			
		def generate_external_scoring(metrics):
			scoring = list()			

			for mclass, runs in metrics.items():
				if mclass in {"repeat", "junction"}:
					continue
				suffixes = NF_SUFFIXES.get(mclass, ["nF1", "jF1", "eF1", "aF1"])
				for runid, run in runs.items():
					for suffix in suffixes:
						comment = "#" if suffix in ["nF1", "jF1", "eF1"] else ""
						multiplier = run[0]["multiplier"]
						if suffix == "tpm":
							expression = "{{rescaling: max, multiplier: {}}}".format(multiplier)
							suffix = ""
						else:
							expression = "{{rescaling: max, use_raw: true, multiplier: {}}}".format(multiplier)
							suffix = "_" + suffix
						scoring.append(
							"  {}external.{}{}: {}".format(comment, runid, suffix, expression)
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


def parse_busco_levels(levels):
	if not levels:
		return True, True, True # while developing
		return True, False, False
	levels = set(l.lower() for l in levels.split(","))
	invalid_levels = levels.difference(BUSCO_LEVELS)
	if invalid_levels:
		raise ValueError("Invalid busco levels specified with --busco-level option ({}). Valid levels are {}.".format(invalid_levels, BUSCO_LEVELS))
	switch_off = {"off", "none"}.intersection(levels)
	unique_levels = set(level[0] for level in levels.difference({"off", "none"}))
	do_proteins = unique_levels.intersection({"a", "p"}) and not switch_off
	do_transcripts = unique_levels.intersection({"a", "t"}) and not switch_off
	do_genome = unique_levels.intersection({"a", "g"}) and not switch_off

	return tuple(map(bool, (do_proteins, do_transcripts, do_genome)))



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
	
	print(args)
	mikado_config_file = os.path.join(args.outdir, args.prefix + ".mikado_config.yaml")
	gmc_config = yaml.load(open(args.config_file), Loader=yaml.SafeLoader)


	cmd = MIKADO_CONFIGURE_CMD.format(
		cmd=gmc_config["program_calls"]["mikado"].format(container=args.mikado_container, program="configure"),
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

	precomputed_busco_genome = dict()
	if args.busco_genome_run is not None:
		try:
			precomputed_busco_genome["summary"] = glob.glob(os.path.join(args.busco_genome_run, "*short_summary*.txt"))[0]
		except:
			raise ValueError("No valid busco genome summary file in {}".format(args.busco_genome_run))
		try:
			precomputed_busco_genome["full_table"] = glob.glob(os.path.join(args.busco_genome_run, "full_table.tsv"))[0]
		except:
			raise ValueError("No valid busco genome full table in {}".format(args.busco_genome_run))
		try:
			precomputed_busco_genome["missing_busco_list"] = glob.glob(os.path.join(args.busco_genome_run, "missing_busco_list.tsv"))[0]
		except:
			raise ValueError("No valid missing busco list in {}".format(args.busco_genome_run))
		
		

	run_config = {
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
		"genus_identifier": args.genus_identifier,
		"busco_analyses": dict(zip(("proteins", "transcriptome", "genome"), parse_busco_levels(args.busco_level))),
		"hpc_config": args.hpc_config
	}

	if any(run_config["busco_analyses"].values()) and (args.busco_lineage is None or not os.path.exists(args.busco_lineage)):
		raise ValueError("BUSCO analysis requested (P:{}, T:{}, G:{}) but no valid lineage specified ({}).".format(
			*run_config["busco_analyses"].values(), 
			args.busco_lineage
		))

	if precomputed_busco_genome and run_config["busco_analyses"]["genome"]:
		raise ValueError("BUSCO genome was selected ({}) together with --busco_genome_run ({}).".format(args.busco_level, args.busco_genome_run))

	run_config["busco_analyses"]["lineage"] = args.busco_lineage	
	run_config["busco_analyses"]["precomputed_genome"] = precomputed_busco_genome
	
	run_data = { 
		"data": {
			"expression-runs": smm.getMetricsData("expression"), 
			"transcript-runs": smm.getMetricsData("aln_tran"), 
			"protein-runs": smm.getMetricsData("aln_prot"), 
			"protein-seqs": smm.getMetricsData("seq_prot"),
			"junction-data": smm.getMetricsData("junction"),	
			"repeat-data": smm.getMetricsData("repeat")
		},
		"transcript_models": {row[1]:row[0] for row in csv.reader(open(args.list_file), delimiter="\t")}
	}


	
	with open(os.path.join(args.outdir, args.prefix + ".run_config.yaml"), "wt") as run_config_out:
		yaml.dump(run_config, run_config_out, default_flow_style=False)
		yaml.dump(run_data, run_config_out, default_flow_style=False)
		yaml.dump(gmc_config, run_config_out, default_flow_style=False)

	pass

