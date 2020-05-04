import csv
from collections import Counter

def generate_final_table(seq_table, bt_conf_table, stats_table, final_table, summary):
	colheaders = ["Confidence", "Biotype", "InFrameStop", "Partialness"]
	pt_cats = {".": "complete", "5_3": "fragment", "3": "3prime_partial", "5": "5prime_partial"}

	if_pt = dict((row[7], row[12:14]) for row in csv.reader(open(seq_table), delimiter="\t"))
	genes, transcripts = dict(), dict()
	for row in csv.reader(open(bt_conf_table), delimiter="\t"):
		if not row[0].startswith("#"):
			if row[2].lower() in {"gene", "mrna", "ncrna"}:
				attrib = dict(item.split("=") for item in row[8].strip().split(";"))
				if row[2].lower() == "gene":
					genes[attrib["ID"]] = (attrib["biotype"], attrib["confidence"])
				else:
					transcripts[attrib["ID"]] = (genes.get(attrib["Parent"], (None, None))[0], attrib["confidence"])

	r = csv.reader(open(stats_table), delimiter="\t")
	head = ["#{}.{}".format(c, col) for c, col in enumerate(next(r), start=1)]
	head.extend("#{}.{}".format(c, col) for c, col in enumerate(colheaders, start=len(head)+1))
	with open(final_table, "w") as tbl_out:
		print(*head, sep="\t", flush=True, file=tbl_out)
		for row in r:
			row.extend(transcripts.get(row[0], [".", "."])[::-1])
			row.extend(if_pt.get(row[0], [".", "."]))
			row[-1] = pt_cats.get(row[-1], "NA")
			print(*row, sep="\t", flush=True, file=tbl_out)

	tcounts, gcounts = Counter(transcripts.values()), Counter(genes.values())
	with open(summary, "w") as summary_out:
		print("Biotype", "Confidence", "Gene", "Transcript", sep="\t", file=summary_out)
		categories = set(tcounts).union(set(gcounts))
		for cat in sorted(categories, key=lambda x:gcounts[x], reverse=True):
			print(*cat, gcounts[cat], tcounts[cat], sep="\t", file=summary_out)
		print("Total", "", sum(gcounts.values()), sum(tcounts.values()), sep="\t", file=summary_out)
