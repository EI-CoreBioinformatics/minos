import sys
import csv
import argparse

CATEGORIES = ["Complete_1", "Complete_2", "Complete_3", "Complete_4+", "Complete", "Duplicated", "Fragmented", "Missing", "Total"]
MAX_COPY_NUMBER = 4

def read_full_table(table_file, tx2gene=None, is_pick=False):
	counts = {cat: 0 for cat in CATEGORIES}
	complete = dict()
	for row in csv.reader(open(table_file), delimiter="\t"):
		if not row or row[0].startswith("#"):
			continue
		if row[1] in {"Missing", "Fragmented"}:
			counts[row[1]] += 1
		else:
			# for duplicated/complete we build lists of (gene) ids to determine the copy number
			if tx2gene is not None:
				# if we have a transcript-gene map (for transcript sets that were processed with mikado prepare)
				gene = tx2gene[row[2]]
			elif is_pick:
				# if this is a final run (on mikado pick output), we assume that transcript_id = gene_id.transcript_number
				gene = row[2][:row[2].rfind(".")]
			else:
				# for genome runs we just add the scaffold/contig
				gene = row[2]
				
			complete.setdefault(row[0], list()).append(tx2gene[row[2]] if tx2gene is not None else row[2])
	for cat, genes in complete.items():
		if tx2gene is not None or is_pick:
			# for non genome runs, check the unique gene ids
			n = len(set(genes))
		else:
			# for genome runs, just count
			n = len(genes)
		cat = "Complete_{}{}".format(min(n, MAX_COPY_NUMBER), "+" if n >= MAX_COPY_NUMBER else "") 
		counts[cat] += 1
		# complete buscos with exactly one unique gene id (or a count of 1) will count towards "Complete", 2+ copies count towards "Duplicated"
		cat = "Complete" if n == 1 else "Duplicated"
		counts[cat] += 1
	# the total in counts.values() contains 2x (duplicate + complete)
	counts["Total"] = sum(v for k, v in counts.items() if not k.startswith("Complete_"))
	return counts

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
