import csv
import re

def generate_tx2gene_maps(_in, _out, runid):
	with open(_out, "w") as tx2gene_out:
		for row in csv.reader(open(_in), delimiter="\t"):
			if row and not row[0].startswith("#"):
				if row[2].lower() in {"mrna", "ncrna"}:
					attr = dict((item.group(1).strip(), item.group(2).strip()) for item in re.finditer("([^;]+)\s*=\s*([^;]+);?", row[8]))
					print("{}_{}".format(runid, attr["ID"]), attr["Parent"], file=tx2gene_out, flush=True, sep="\t")
