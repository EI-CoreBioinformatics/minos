import csv
import sys
import yaml
import argparse
import os
import pathlib
import subprocess

from collections import OrderedDict

from . import __version__

def parseListFile(fn):
	d = OrderedDict()
	for row in csv.reader(open(fn), delimiter="\t"):
		d[row[1]] = row[0], bool(row[2]), int(row[3]), bool(row[4])
	return d

def createScoringFile(fn, hints, fo):
	# gather external hints 
	coding = list()
	for k in hints:
		if k.endswith("_coding"):
			coding.append(k.strip("_coding"))
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


def run_configure(args):
	listFile = parseListFile(args.list_file)

	pathlib.Path(args.outdir).mkdir(exist_ok=True, parents=True)
	scoringFile = os.path.join(args.outdir, args.prefix + ".scoring.yaml")

	createScoringFile(args.scoring_template, listFile, scoringFile)

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

	pass

