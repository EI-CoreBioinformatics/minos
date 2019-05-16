import sys
import argparse
import csv
import collections




def main():

	ap = argparse.ArgumentParser()
	ap.add_argument("gff", type=str)
	ap.add_argument("--source", type=str, default="Mikado_annotation_run2")
	args = ap.parse_args()

	for row in csv.reader(open(args.gff), delimiter="\t"):
		if row[0]:
			if not row[0].startswith("#"):
				feature = row[2]
				row[1] = args.source
				if feature == "sublocus":
					row[2] = "gene"
					row[8] = row[8].replace("Parent=", "superlocus=")
					print("###")
				elif feature == "ncRNA_gene":
					row[2] == "gene"
				elif feature in {"mRNA", "ncRNA", "transcript"}:
					row[2] = "mRNA"
					attrib = dict(p.split("=") for p in row[8].split(";"))
					tid = attrib.get("ID", "NA")
					try:
						alias = attrib["alias"]
					except:
						print(
							"# alias attribute not defined for transcript {}, checking target instead ... ".format(tid), 
							end="", 
							file=sys.stderr
						)
						try:
							alias = attrib["target"]
							print("# Got target, using it as Name attribute.", file=sys.stderr)
						except:
							print(
								"# target attribute not defined for transcript {}, checking prev_parent instead ... ".format(tid),
								end="",
								file=sys.stderr)
							try:
								alias = attrib["prev_parent"]
								print("# Got prev_parent, using it as Name attribute", file=sys.stderr)
							except:
								print(
									"# prev_parent attribute not defined for transcript {}. Falling back to ID ... ".format(tid),
									end="",
									file=sys.stderr
								)
								alias = tid
					preserve_attribs = ("ID", "Parent", "Name", "Note")
					new_attrib = collections.OrderedDict((k, attrib.get(k)) for k in preserve_attribs)
					new_attrib["Note"] = [tid]
					new_attrib["Note"].extend("{}:{}".format(k, attrib[k]) for k in sorted(attrib) if k not in preserve_attribs)
					row[8] = ";".join("{}={}".format(k, new_attrib[k]) for k in new_attrib)
									
			print(*row, sep="\t")
					 


			
			
	



	pass


if __name__ == "__main__": 
	main()
