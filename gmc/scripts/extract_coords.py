import csv
import re

def extract_coords(_in, _out, filetype="gtf"):
	regex = re.compile('transcript_id "?([^";]+)"?;' if filetype == "gtf" else 'ID\s?=\s?"?([^";]+)"?;')
	with open(_out, "w") as coords_out:
		for row in csv.reader(open(_in), delimiter="\t"):
			if row and not row[0].startswith("#") and "rna" in row[2].lower():
				tid = regex.search(row[8]).group(1)
				print(tid, "{seq}:{start}..{end}".format(seq=row[0], start=row[3], end=row[4]), sep="\t", flush=True, file=coords_out)
