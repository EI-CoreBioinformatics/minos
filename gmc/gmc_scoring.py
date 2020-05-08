import csv
import sys
import os

from collections import OrderedDict

EXTERNAL_METRICS_HEADERS = ["metric_name_prefix", "metric_class", "multiplier", "not_fragmentary_min_value", "file_path"]

NF_SUFFIXES = {
	"seq_prot": ["qCov", "tCov"],
	"blastdb_prot": ["qCov", "tCov"],
	"repeat": ["cov"],
	"expression": ["tpm"],
	"busco": ["busco"]
}

NON_NF = {"junction", "expression", "repeat", "busco"}

METRIC_KEY_LABELS = {
	"expression": "expression-runs",
	"aln_tran": "transcript-runs",
	"aln_prot": "protein-runs",
	"seq_prot": "protein-seqs",
	"junction": "junction-data",
	"repeat": "repeat-data"
}




class ScoringMetricsManager(object):
	STRANDINFO = {"_xx": "unstranded", "_rf": "rf-stranded", "_fr": "fr-stranded"}
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
					if metric_name[-3:] not in ScoringMetricsManager.STRANDINFO:
						raise ValueError("Expression metric {} does not have strandedness information. Please add a suffix to indicate strandedness ({})".format(
							metric_name, ScoringMetricsManager.STRANDINFO
						))
				else:
					print("Found expression metric: {} but --use-tpm-for-picking was not set. Ignoring.".format(metric_name))
					continue

			self.metrics.setdefault(metric_class, OrderedDict()).setdefault(metric_name, list()).append(data)
	
		if len(self.metrics.get("junction", OrderedDict())) > 1:
			raise ValueError("More than one junction file supplied. " + self.metrics.get("junction", list()))


	def __init__(self, metrics_file, scoring_template_file, outdir, prefix, use_tpm=False):
		self.__importMetricsData(metrics_file, use_tpm=use_tpm)

	def getMetricsData(self, metric):
		mdata = dict()
		for k in self.metrics.get(metric, OrderedDict()):
			mdata.setdefault(k, list()).extend(item["files"] for item in self.metrics[metric][k])
			pass
		return mdata
	
	def generateScoringFile(self, scoring_template, outfile):

		def generate_nf_expression(metrics):
			ext_metrics = list()
			for mclass, runs in self.metrics.items():
				if mclass not in NON_NF:
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
				if mclass not in NON_NF:
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
					print("  external.busco: {rescaling: max, use_raw: true, multiplier: 1}", file=_out)
					# print("  # all boolean metrics values from here below", file=_out)
					# print("  # external.EI_tpm_05: {rescaling: max, use_raw: true, multiplier: 10}", file=_out)
					# print("  # external.EI_tpm_1: {rescaling: max, use_raw: true, multiplier: 10}", file=_out)
					# print("  # external.SRA_tpm_05: {rescaling: max, use_raw: true, multiplier: 10}", file=_out)
					# print("  # external.SRA_tpm_1: {rescaling: max, use_raw: true, multiplier: 10}", file=_out)
					# print("  external.fln: {rescaling: max, use_raw: true, multiplier: 5}", file=_out)
				else:
					print(line, end="", file=_out)				

	def get_scoring_data(self):
		return {label: self.getMetricsData(key) for key, label in METRIC_KEY_LABELS.items()}
