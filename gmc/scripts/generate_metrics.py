import sys
import os
import csv
import argparse
import yaml
import collections


def read_metrics_info(f):
	metrics_info = collections.OrderedDict()
	for row in csv.reader(open(f), delimiter="\t"):
		metrics_type = "mikado" if "mikado" in row[0] else row[0]
		metrics_info.setdefault(metrics_type, collections.OrderedDict())[row[1]] = row[2]

	return metrics_info

def generate_metrics_header(metrics_info):
	header = ["tid"]
	for mcid in MCLASS_INFO:
		for mid in metrics_info.get(mcid, collections.OrderedDict()):
			header.extend([mid + (("_" + m) if m else m) for m in MCLASS_INFO[mcid]["metrics"]])

	return header

def _round_perc_frac(x):
	return '%.4f' % (x / 100)

def _round_frac(x):
	return '%.4f' % x

def clean_na(x):
	try:
		return float(x)
	except:
		return 0.0


def read_mikado_refmap(f, seen=set()):
	model_info = {"nmetrics": None}
	for i, row in enumerate(csv.reader(open(f), delimiter="\t")):
		if i > 0 and not row[0].startswith("#"):
			tid = row[0]
			if seen and tid not in seen:
				raise ValueError("Error: This is not the first processed mikado run ({}) but transcript {} occurs for the first time. Please check your mikado inputs.".format(f, tid))

			nF1, jF1, eF1 = map(clean_na, (row[6], row[9], row[12]))
			aF1 = sum((nF1, jF1, eF1)) / 3
			model_info[tid] = collections.OrderedDict(zip(("nF1", "jF1", "eF1", "aF1"), map(_round_perc_frac, (nF1, jF1, eF1, aF1))))

	if seen:
		delta = seen.difference(model_info)
		if delta:
			raise ValueError("Error: Mikado run ({}) has less transcripts than previously processed run ({} missing transcript(s): {}). Please check your mikado inputs.".format(f, len(delta), ",".join(sorted(delta)[:3]) + ", ..."))

	return model_info

def read_blast(f, seen=set()):
	model_info = {"nmetrics": 2}
	for i, row in enumerate(csv.reader(open(f), delimiter="\t")):
		if not row[0].startswith("#"):
			tid = row[0]
			if seen and tid not in seen:
				raise ValueError("Error: This is a blast run ({}) but transcript {} occurs for the first time. Please check your mikado/blast inputs.".format(f, tid))

			qCov, tCov = map(clean_na, (row[17], row[18]))
			model_info[tid] = collections.OrderedDict(zip(("qCov", "tCov"), map(_round_perc_frac, (qCov, tCov))))
			
	return model_info

def read_cpc(f, seen=set()):
	model_info = {"nmetrics": None}
	for i, row in enumerate(csv.reader(open(f), delimiter="\t")):
		if i > 0 and not row[0].startswith("#"):
			tid = row[0]
			if seen and tid not in seen:
				raise ValueError("Error: This is a cpc run ({}) but transcript {} occurs for the first time. Please check your mikado/cpc inputs.".format(f, tid))
			
			model_info[tid] = collections.OrderedDict({"cpc": _round_frac(clean_na(row[6]))})

	return model_info

def read_kallisto(f, seen=set()):
	model_info = {"nmetrics": None}
	for i, row in enumerate(csv.reader(open(f), delimiter="\t")):
		if i > 0 and not row[0].startswith("#"):
			tid = row[0]
			if seen and tid not in seen:
				raise ValueError("Error: This is a kallisto run ({}) but transcript {} occurs for the first time. Please check your mikado/kallisto inputs.".format(f, tid))

			model_info[tid] = collections.OrderedDict({"tpm": _round_frac(clean_na(row[4]))})

	return model_info

def read_rm_repeats(f, seen=set()):
	model_info = {"nmetrics": 1}
	for i, row in enumerate(csv.reader(open(f), delimiter="\t")):
		if row and not row[0].startswith("#"):
			tid = row[0]
			if seen and tid not in seen:
				raise ValueError("Error: This is a repeat modeler run ({}) but transcript {} occurs for the first time. Please check your mikado/rm inputs.".format(f, tid))

			model_info[tid] = collections.OrderedDict({"cov": _round_frac(clean_na(row[3]))})

	return model_info

def read_busco(f, seen=set(), score_complete=1, score_fragmented=0.25):
	model_info = {"nmetrics": 1}
	for i, row in enumerate(csv.reader(open(f), delimiter="\t")):
		if row and not row[0].startswith("#"):
			if row[1] != "Missing":
				tid = row[2]
				if seen and tid not in seen:
					raise ValueError("Error: This is a busco run ({}) but transcript {} occurs for the first time. Please check your mikado/rm inputs.".format(f, tid))

				score = score_complete if row[1] in {"Duplicated", "Complete"} else score_fragmented if row[1] == "Fragmented" else 0
				model_info[tid] = collections.OrderedDict({"busco": _round_frac(clean_na(score))})

	return model_info


MCLASS_INFO = collections.OrderedDict([
	("mikado", {"metrics": ["nF1", "jF1", "eF1", "aF1"], "parser": read_mikado_refmap}),
	("blast", {"metrics": ["qCov", "tCov"], "parser": read_blast}),
	("cpc", {"metrics": [""], "parser": read_cpc}),
	("expression", {"metrics": [""], "parser": read_kallisto}),
	("repeat", {"metrics": ["cov"], "parser": read_rm_repeats}),
	("busco", {"metrics": ["busco"], "parser": read_busco})	
])                                                                                      		


def main():

	ap = argparse.ArgumentParser()
	ap.add_argument("metrics_info", type=str)
	args = ap.parse_args()

	if not os.path.exists(args.metrics_info):
		raise ValueError("Error: Could not find metrics info at " + args.metric_info)

	metrics_info = read_metrics_info(args.metrics_info)
	metrics_header = generate_metrics_header(metrics_info)

	print(*metrics_header, sep="\t")

	model_info = collections.OrderedDict()
	first_run = None
	seen = set()
	for mclass in MCLASS_INFO:
		for run in metrics_info.get(mclass, collections.OrderedDict()):
			parser = MCLASS_INFO[mclass]["parser"]
			if parser is not None:
				if first_run is None:
					first_run = run
				model_info[run] = parser(metrics_info[mclass][run], seen=seen)
				if not seen:
					seen = set(model_info[run])

	for tid in model_info[first_run]:
		if tid != "nmetrics":
			row = [tid]
			for run in model_info:
				try:
					data = model_info[run][tid].values()
				except KeyError:
					try:
						data = ["0.0000" for i in range(model_info[run]["nmetrics"])]
					except ValueError:
						raise ValueError("Error: Missing {} values for transcript {}".format(run, tid))
				row.extend(data)
			print(*row, sep="\t")


if __name__ == "__main__":
	main()
