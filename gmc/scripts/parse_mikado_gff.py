import sys
import argparse
import csv
import collections


def get_alias(tid, attrib):

	err_msg = "# {} attribute not defined for transcript {}. {}"
	resolution = "Checking {} instead ... "
	success_msg = "# Got {}, using it as Name attribute."

	try:
		alias = attrib["alias"]
	except:
		print(err_msg.format("alias", tid, resolution.format("target")), end="", file=sys.stderr)
		try:
			alias = attrib["target"]
			print(success_msg.format("target"), file=sys.stderr)
		except:
			print(err_msg.format("target", tid, resolution.format("prev_parent")), end="", file=sys.stderr)
			try:
				alias = attrib["prev_parent"]
				print(success_msg.format("prev_parent"), file=sys.stderr)
			except:
				print(err_msg.format("prev_parent", tid, "Falling back to ID ... "), end="", file=sys.stderr)
				alias = tid

	return alias


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
				if feature == "superlocus":
					continue
				if feature == "sublocus":
					row[2] = "gene"
					row[8] = row[8].replace("Parent=", "superlocus=")
					print("###")
				elif feature == "ncRNA_gene":
					row[2] == "gene"
				elif feature in {"mRNA", "ncRNA", "transcript"}:
					row[2] = "mRNA"
					attrib = dict(p.split("=") for p in row[8].split(";"))
					tid = attrib.get("ID", None)
					if tid is None:
						raise ValueError("Could not find transcript id.\n{}\n".format(*row, sep="\t"))

					alias = get_alias(tid, attrib)
					preserve_attribs = ("ID", "Parent", "Name", "alias", "Note")
					new_attrib = {
						"ID": tid,
						"Parent": attrib["Parent"],
						"Name": alias,
						"Note": [tid]
					}

					# new_attrib = collections.OrderedDict((k, attrib.get(k)) for k in preserve_attribs)
					# new_attrib["alias"] = alias
					# new_attrib["Note"] = [tid]
					new_attrib["Note"].extend("{}:{}".format(k, v) for k, v in attrib.items() if k not in preserve_attribs and k.lower() != "note")
					new_attrib["Note"].extend(attrib.get("note", "").split("|")[1:])
					new_attrib["Note"] = ",".join(new_attrib["Note"])
					row[8] = ";".join("{}={}".format(k, v) for k, v in new_attrib.items())
				
					
			print(*row, sep="\t")
					 


			
			
	



	pass


if __name__ == "__main__": 
	main()
