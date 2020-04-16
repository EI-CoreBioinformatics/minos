
def split_fasta(fastafile, split_files):
	out = None
	for line in open(fastafile):
		if line.startswith(">"):
			matches = list({runid for runid in split_files if line[1:].startswith(runid)})
			if not matches:
				print("No matching output file for sequence " + line[1:].strip())
				out = None
				continue
			if len(matches) > 1:
				matches = sorted(matches, key=lambda x:len(x), reverse=True)
			out = split_files[matches[0]]
		if out is not None:
			print(line, end="", file=out, flush=True)

	for f in split_files.values():
		f.close()
