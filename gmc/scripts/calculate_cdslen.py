def calculate_cdslen(_in, _out, min_cds_length):
	def read_fasta(f):
		sid, seq = None, list()
		for line in open(f):
			line = line.strip()
			if line.startswith(">"):
				if seq:
					yield sid, "".join(seq)
					seq = list()
				sid = line[1:].split()[0]
			else:
				seq.append(line)
		if seq:
			yield sid, "".join(seq)

	with open(_out, "w") as fout:
		for sid, seq in read_fasta(_in):
			seq = seq.upper().replace("U", "T")
			has_stop = any(seq.endswith(stop_codon) for stop_codon in ("TAA", "TGA", "TAG"))
			discard = (has_stop and len(seq) < min_cds_length) or (not has_stop and len(seq) < min_cds_length - 3)
			print(sid, int(discard), sep="\t", flush=True, file=fout)
