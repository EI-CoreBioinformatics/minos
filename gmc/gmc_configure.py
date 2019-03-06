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
	def __importMetricsData(self, fn):
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
				if row[0][-3:] not in STRANDINFO:
					raise ValueError("ERROR: expression metric does not have strandedness information " + row[0])				

			self.metrics.setdefault(row[1], OrderedDict()).setdefault(row[0], list()).append(data)


		for k in self.metrics:
			print(k)
			for j in self.metrics[k]:
				print(j, self.metrics[k][j], sep="\n")

	def __init__(self, metrics_file, scoring_template_file, outdir, prefix):
		self.__importMetricsData(metrics_file)

	"""
EI_kallisto_rf
[
{'multiplier': '0', 'min_value': '0', 'files': ['/ei/workarea/group-ga/Projects/CB-GENANNO-430_Annotation_of_Melia_azedarach_and_Quillaja_saponaria/Analysis/Quillaja_saponaria/trim_galore-0.5.0_gk/Expanding_leaf/Expanding_leaf_R1_val_1.fq.gz', '/ei/workarea/group-ga/Projects/CB-GENANNO-430_Annotation_of_Melia_azedarach_and_Quillaja_saponaria/Analysis/Quillaja_saponaria/trim_galore-0.5.0_gk/Expanding_leaf/Expanding_leaf_R2_val_2.fq.gz']}, 
{'multiplier': '0', 'min_value': '0', 'files': ['/ei/workarea/group-ga/Projects/CB-GENANNO-430_Annotation_of_Melia_azedarach_and_Quillaja_saponaria/Analysis/Quillaja_saponaria/trim_galore-0.5.0_gk/Green_stem/Green_stem_R1_val_1.fq.gz', '/ei/workarea/group-ga/Projects/CB-GENANNO-430_Annotation_of_Melia_azedarach_and_Quillaja_saponaria/Analysis/Quillaja_saponaria/trim_galore-0.5.0_gk/Green_stem/Green_stem_R2_val_2.fq.gz']}
]
	"""

	def getMetricsData(self, metric):
		mdata = dict()
		for k in self.metrics.get(metric, OrderedDict()):
			mdata.setdefault(k, list()).extend(item["files"] for item in self.metrics[metric][k])
			pass
		return mdata
	
	def generateScoringFile(self, scoring_template_file, outfile):
		metrics_ids = ["external.{}_aF1".format(k) for k in self.metrics]
		
		with open(scoring_template_file) as _in, open(outfile, "wt") as _out:
			print(_in.read(), end="", file=_out)
			"""
			for line in _in:
				print(line, end="", file=_out)
				if line.strip().startswith("not_fragmentary:"):
					break

				expression = "[((exon_num.multi and (combined_cds_length.multi or {0}))" + \
					", or, " + \
					"(exon_num.mono and (combined_cds_length.mono or {0})))]"
				expression = expression.format("*".join(["external.all_aF1"] + metrics_ids).replace("*", " or "))
				print("  expression: " + expression, file=_out)
				for line in _in:
					if line.strip().startswith("expression:"):
						line = line.replace("expression:", "# expression:")
					print(line, end="", file=_out)
					if line.strip().startswith("external.all_aF1"):
						for m in metrics_ids:
							print(line.replace("external.all_aF1", m), end="", file=_out)
					if line.strip().endswith("external metrics START"):
						break

			for line in _in:
				print(line, end="", file=_out)
			"""


def parseListFile(fn):
	d = OrderedDict()
	for row in csv.reader(open(fn), delimiter="\t"):
		d[row[1]] = row[0], bool(row[2]), int(row[3]), bool(row[4])
	return d

def createScoringFile(fn, hints, fo):
	# gather external hints 
	print("Processing list data:")
	coding = list()
	for k in hints:
		print(k, k.replace("_coding", "") if k.endswith("_coding") else ".", sep="\t")
		if k.endswith("_coding"):
			coding.append(k.replace("_coding", ""))
	metrics = ["external.{}_aF1".format(k) for k in coding]	

	# parse template
	with open(fn) as _in, open(fo, "wt") as _out:
		for line in _in:
			print(line, end="", file=_out)
			if line.strip().startswith("not_fragmentary:"):
				break

		expr = "[((exon_num.multi and (combined_cds_length.multi or {0}))" + \
			", or, " + \
			"(exon_num.mono and (combined_cds_length.mono or {0})))]"
		expr = expr.format("*".join(["external.all_aF1"] + metrics).replace("*", " or "))
		print("  expression: " + expr, file=_out)
		for line in _in:
			if line.strip().startswith("expression:"):
				line = line.replace("expression:", "# expression:")
				
			print(line, end="", file=_out)
			if line.strip().startswith("external.all_aF1"):
				for m in metrics:
					print(line.replace("external.all_aF1", m), end="", file=_out)
			if line.strip().endswith("external metrics START"):
				break


		for m in ["external.all_aF1", "external.mikado_aF1"] + metrics:
			for sfx in ["nF1", "jF1", "eF1", "aF1"]:
				multiplier = 10 if sfx == "aF1" else (5 if not "mikado" in m else 2)
				comment = "# " if not sfx == "aF1" else ""
				print("  " + comment + m.replace("_aF1", "_" + sfx) + ": {{rescaling: max, use_raw: true, multiplier: {}}}".format(multiplier), file=_out)
		for m in metrics:
			for sfx in ["qCov", "tCov"]:
				print("  " + m.replace("_aF1", "_" + sfx) + ": {rescaling: max, use_raw: true, multiplier: 5}", file=_out)

		for line in _in:
			print(line, end="", file=_out)


def parse_external_metrics(fn):
	expression_runs = dict()
	transcript_runs = dict()
	protein_runs = dict()

	for row in csv.reader(open(fn), delimiter="\t", quotechar="\""):
		if row[0].startswith("#"):
			continue
		# metric_name_prefix    metric_class    multiplier  not_fragmentary_min_value   file_path		
		if row[1] == "expression":
			if row[0][-3:] not in {"_xx", "_rf", "_fr"}:
				raise ValueError("ERROR: expression metric does not have strandedness information " + row[0])
			expression_runs.setdefault(row[0], list()).append(row[4].split(","))
		elif row[1] == "aln_tran":
			transcript_runs.setdefault(row[0], list()).append(row[4])
		elif row[1] == "aln_prot": 
			protein_runs.setdefault(row[0], list()).append(row[4])

	return expression_runs, transcript_runs, protein_runs
			
		


def run_configure(args):

	pathlib.Path(args.outdir).mkdir(exist_ok=True, parents=True)
	pathlib.Path(os.path.join(args.outdir, "hpc_logs")).mkdir(exist_ok=True, parents=True)

	scoringFile = os.path.join(args.outdir, args.prefix + ".scoring.yaml")
	smm = ScoringMetricsManager(args.external_metrics, args.scoring_template, args.outdir, args.prefix)
	print("Generating scoring file " + scoringFile + " ...", end="", flush=True)
	smm.generateScoringFile(args.scoring_template, scoringFile)
	print(" done.")


	# parse external metrics here
	#expression_runs, transcript_runs, protein_runs = dict(), dict(), dict()
	#if args.external_metrics:
	#	expression_runs, transcript_runs, protein_runs = parse_external_metrics(args.external_metrics)


	#listFile = parseListFile(args.list_file)
	#createScoringFile(args.scoring_template, listFile, scoringFile)



	#!TODO: 
	# - scan config template for reference
	# - warn if not present
	# - add command line option 
	
	mikado_config_file = os.path.join(args.outdir, args.prefix + ".mikado_config.yaml")

	cmd = "singularity exec {} mikado configure --list {} {} -od {} --scoring {} {}".format(
		args.mikado_container,
		args.list_file,
		("--external " + args.external) if args.external else "",
		args.outdir,
		scoringFile,
		# os.path.abspath("mikado_config.yaml")
		mikado_config_file
	)

	print(cmd)
	out = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)

	#with open(mikado_config_file, "wt") as config_out:
	#	print(out.decode(), sep="\n", file=config_out)


	run_zzz_config = {
		"prefix": args.prefix,
		"outdir": args.outdir,
		"mikado-container": args.mikado_container,
		"mikado-config-file": mikado_config_file,
		"external-metrics": args.external_metrics,
		"expression-runs": smm.getMetricsData("expression"), #expression_runs,
		"transcript-runs": smm.getMetricsData("aln_tran"), #transcript_runs,
		"protein-runs": smm.getMetricsData("aln_prot") #protein_runs
	}
	
	with open(os.path.join(args.outdir, args.prefix + ".run_config.yaml"), "wt") as run_config_out:
		yaml.dump(run_zzz_config, run_config_out, default_flow_style=False)

	pass

