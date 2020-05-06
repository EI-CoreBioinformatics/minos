import sys
import csv
import argparse

from collections import Counter

MAIN_CATEGORIES = ["Complete", "Duplicated", "Fragmented", "Missing", "Total"]

def get_busco_categories(max_copy_number=4):
	return ["Complete_{}{}".format(n, "+" if n == max_copy_number else "") for n in range(1, max_copy_number + 1)] + MAIN_CATEGORIES


def read_full_table(table_file, tx2gene=None, is_pick=False, max_copy_number=4):

	def get_gene_id(seqid, tx2gene=None, is_pick=False):
		if tx2gene is not None:
			# if we have a transcript-gene map (for transcript sets that were processed with mikado prepare)
			return tx2gene[seqid]
		elif is_pick:
			# if this is a final run (on mikado pick output), we assume that transcript_id = gene_id.transcript_number
			return seqid[:seqid.rfind(".")]
		else:
			# for genome runs we just add the scaffold/contig
			return seqid

	counts = Counter({cat: 0 for cat in get_busco_categories(max_copy_number=max_copy_number)})
	complete, missing, fragmented = dict(), set(), dict()
	for row in csv.reader(open(table_file), delimiter="\t"):
		if not row or row[0].startswith("#"):
			continue
		if row[1] in {"Missing", "Fragmented"}:
			if row[1] == "Missing":
				missing.add(row[0])
			else:
				fragmented.setdefault(row[0], list()).append(row[2])
			counts[row[1]] += 1
		else:
			# for duplicated/complete we build lists of (gene) ids to determine the copy number
			complete.setdefault(row[0], list()).append((row[2], get_gene_id(row[2], tx2gene=tx2gene, is_pick=is_pick)))

	for cat, genes in complete.items():
		if tx2gene is not None or is_pick:
			# for non genome runs, check the unique gene ids
			n = len(set(gid for tid, gid in genes))
		else:
			# for genome runs, just count
			n = len(genes)
		cat = "Complete_{}{}".format(min(n, max_copy_number), "+" if n >= max_copy_number else "") 
		counts[cat] += 1
		# all complete buscos will count towards "Complete", 2+ copies count towards "Duplicated"
		if n > 1:
			counts["Duplicated"] += 1

	counts["Complete"] = sum(v for k, v in counts.items() if k.startswith("Complete_"))
	counts["Total"] = counts["Complete"] + counts["Fragmented"] + counts["Missing"]
	return counts, complete, missing, fragmented

def read_tx2gene(f):
	return dict(row[:2] for row in csv.reader(open(f), delimiter="\t"))

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("table_file", type=str)
	ap.add_argument("tx2gene", type=str)
	args = ap.parse_args()

	tx2gene = read_tx2gene(args.tx2gene)
	busco_table = read_full_table(args.table_file, tx2gene)	
	print(busco_table)


if __name__ == "__main__":
	main()
