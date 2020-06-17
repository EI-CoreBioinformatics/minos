import os
import glob
import sys
import csv
import re

from minos.minos_configure import ExternalMetrics

def generate_metrics_info(metrics_path, _out, busco_data, use_diamond):

	# this block is unstable against tempering with the structure of the generate_metrics dir
	# might work better as state-machine
	with open(_out, "wt") as metrics_info:
		walk = os.walk(metrics_path, followlinks=True)
		next(walk)
		rows = list()

		while walk:
			try:
				cwd, dirs, files = next(walk)
			except StopIteration:
				break
			cwd_base = os.path.basename(cwd)
			if cwd_base == "CPC-2.0_beta":
				mclass, mid, path = "cpc", "cpc", os.path.join(cwd, files[0])
				rows.append((ExternalMetrics.CPC_CODING_POTENTIAL, mclass, mid, path))
			elif cwd_base in {"proteins", "transcripts"}:
				mclass = "mikado.{}".format(cwd_base[:-1])
				for mid in dirs:
					cwd, _, files = next(walk)
					path = glob.glob(os.path.join(cwd, "*.refmap"))[0]
					rows.append((ExternalMetrics.MIKADO_TRANSCRIPTS_OR_PROTEINS, mclass, mid, path))
			elif cwd_base in {"blastp", "blastx"}:
				mclass = "blast"
				for mid in dirs:
					cwd, dbdirs, files = next(walk)
					pattern = (r"\DDIAMOND\D*\Dtophit" if use_diamond else r"\DBLAST\D*\Dtophit")
					for f in filter(lambda x: re.search(pattern, x), files):
						path = os.path.join(cwd, f)
						rows.append((ExternalMetrics.PROTEIN_BLAST_TOPHITS, mclass, mid, path))
					# skip the blast databases (if existing)
					for dbdir in dbdirs:
						_ = next(walk)
			elif cwd_base == "kallisto":
				mclass = "expression"
				for mid in dirs:
					cwd, _, _ = next(walk)
					path = os.path.join(cwd, "abundance.tsv")
					rows.append((ExternalMetrics.KALLISTO_TPM_EXPRESSION, mclass, mid, path))
			elif cwd_base == "repeats":
				mclass = "repeat"
				for path in glob.glob(os.path.join(cwd, "*.no_strand.exon.gff.cbed.parsed.txt")):
					mid = os.path.basename(path).split(".")[0]
					rows.append((ExternalMetrics.REPEAT_ANNOTATION, mclass, mid, path))
			elif cwd_base == "busco_proteins":
				mclass = "busco"
				mid = busco_data
				path = os.path.join(metrics_path, "busco_proteins", "busco_proteins.tsv")
				rows.append((ExternalMetrics.BUSCO_PROTEINS, mclass, mid, path))
		"""
		for mid in config["data"].get("repeat-data", dict()):
			mclass = "repeat"
			path = config["data"]["repeat-data"][mid][0][0]
			rows.append((ExternalMetrics.REPEAT_ANNOTATION, mclass, mid, path))
		"""

		for _, mclass, mid, path in sorted(rows, key=lambda x:x[0].value):
			print(mclass, mid, path, sep="\t", file=sys.stderr)
			print(mclass, mid, os.path.abspath(path), sep="\t", file=metrics_info)
