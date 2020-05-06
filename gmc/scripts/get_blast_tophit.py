import csv

def get_blast_tophit(_in, _out, pident_threshold, qcov_threshold):
	def calc_coverage_perc(start, end, length):
		alen = (end - start + 1) if start < end else (start - end + 1)
		return round(alen/length * 100, 2)
	with open(_out, "w") as blast_out:
		seen = set()
		for row in csv.reader(open(_in), delimiter="\t"):
			if row and not row[0].startswith("#") and not row[0] in seen:
				if len(row) < 17:
					raise ValueError("Misformatted blastp line (ncols={ncols}) in {blastp_input}. Expected is outfmt 6 (17 cols).".format(
						ncols=len(row),
						blastp_input=_in)
					)
				try:
					qstart, qend, sstart, send, qlen, slen = map(int, row[3:9])
				except:
					raise ValueError("Misformatted blastp line: could not parse integer values from columns 2-7 ({})".format(",".join(row[3:9])))
				qcov, scov = calc_coverage_perc(qstart, qend, qlen), calc_coverage_perc(sstart, send, slen)
				if float(row[2]) >= pident_threshold and qcov >= qcov_threshold:
					print(*row, "{:.2f}".format(qcov), "{:.2f}".format(scov), sep="\t", file=blast_out)
					seen.add(row[0])
