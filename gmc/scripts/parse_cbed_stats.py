import csv
import sys

def parse_cbed(instream, print_header=False, outstream=sys.stdout):
	def extract_minmax(coords):
		_min, _max = coords[0]
		for start, end in coords[1:]:
			_min, _max = min(_min, start), max(_max, end)
		return _min, _max
	def write_record(scaffold, parent, trans_bps, covered_bps, scaffold_coords, out):
		start, end = extract_minmax(scaffold_coords)
		coverage = (covered_bps / trans_bps) if trans_bps else None
		print(cur_parent, trans_bps, covered_bps, "{:.2f}".format(coverage) if coverage is not None else "NA", "{}:{}..{}".format(scaffold, start, end), sep="\t", file=out)
	if print_header:
		print("#ID", "#Total_bps", "#bps_covered", "#%bps_covered", "#scaffold_BrowserView", sep="\t", file=outstream)
	cur_parent = None
	scaffold_coords = list()
	covered_bps, trans_bps = 0, 0
	scaffold = str()
	for row in csv.reader(instream, delimiter="\t"):
		if row and not row[0].startswith("#"):
			attr = dict(item.split("=") for item in row[8].strip(" ;").split(";"))
			if attr["Parent"] != cur_parent:
				if cur_parent is not None:
					write_record(scaffold, cur_parent, trans_bps, covered_bps, scaffold_coords, outstream)
					scaffold_coords.clear()
					pass
				cur_parent = attr["Parent"]
				covered_bps, trans_bps = 0, 0
				scaffold = row[0]

			covered_bps += int(row[10])
			trans_bps += int(row[11])
			start, end = int(row[3]), int(row[4])
			if end < start:
				start, end = end, start
			scaffold_coords.append((start, end))

	write_record(scaffold, cur_parent, trans_bps, covered_bps, scaffold_coords, outstream)
