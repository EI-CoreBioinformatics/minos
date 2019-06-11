import os
import argparse
from collections import Counter

from Bio import SeqIO


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--outdir", "-o", default=".", type=str)
	ap.add_argument("--prefix", "-p", default="", type=str)


	ap.add_argument("fasta")



	args = ap.parse_args()
	
	with open(os.path.join(args.outdir, args.prefix + ".protein_status.tsv"), "w") as p_out:

		print("#1.protein_id", "#2.status", sep="\t", file=p_out)
		c = Counter()
		for record in SeqIO.parse(args.fasta, "fasta"):
			has_start = record.seq[0].upper() == "M"
			has_stop = record.seq[-1] == "."
			n_stops = record.seq.count(".")
	
			actual_inframe_stops = n_stops - 1 if has_stop else n_stops
			has_inframe_stop = actual_inframe_stops > 0
	
			status = None
			if has_start and has_stop and not has_inframe_stop:
				status = "complete"
			elif has_start and not has_stop and not has_inframe_stop:
				status = "3prime_partial"
			elif not has_start and has_stop and not has_inframe_stop:
				status = "5prime_partial"
			elif has_inframe_stop:
				status = "in-frame_stop"
			elif not has_start and not has_stop and not has_inframe_stop:
				status = "fragment"
			else:
				status = "unknown"
	
			if status is None:
				raise ValueError("Error: Status cannot be determined for sequence " + record.id)

			c[status] += 1
	
			print(record.id, status, sep="\t", file=p_out)


	with open(os.path.join(args.outdir, args.prefix + ".protein_status.summary"), "w") as p_summary_out:
		for k, v in sorted(c.items()):
			print(k, v, sep="\t", file=p_summary_out)

					


if __name__ == "__main__":
	main()
