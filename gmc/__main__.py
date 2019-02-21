import csv
import sys
import yaml
import argparse
import os
import pathlib
import subprocess

from collections import OrderedDict

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



			
		

		



def main():	

	ap = argparse.ArgumentParser()
	ap.add_argument("list_file", type=str)
	ap.add_argument("scoring_template", type=str)

	ap.add_argument("--external", type=str, default="")
	ap.add_argument("--outdir", "-o", type=str, default="gmc_run")
	ap.add_argument("--mikado-container", type=str, default="/ei/software/testing/gmc/dev/x86_64/mikado.simg")

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
	

	listFile = parseListFile(args.list_file)

	pathlib.Path(args.outdir).mkdir(exist_ok=True, parents=True)
	scoringFile = os.path.join(args.outdir, "scoring_test.yaml")

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

	print(out.decode(), sep="\n")

	


	pass

if __name__ == "__main__": 
	main()
