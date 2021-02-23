def calculate_cdslen(_in_cds, _in_cdna, _out, min_cds_length):
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
		cdna_set = set(sid for sid, _ in read_fasta(_in_cdna))
		cds_set = set()
		for sid, seq in read_fasta(_in_cds):
			cds_set.add(sid)
			seq = seq.upper().replace("U", "T")
			has_stop = any(seq.endswith(stop_codon) for stop_codon in ("TAA", "TGA", "TAG"))
			discard = (has_stop and len(seq) < min_cds_length) or (not has_stop and len(seq) < min_cds_length - 3)
			print(sid, int(discard), sep="\t", flush=True, file=fout)

		# get the ncRNA models and add them to CDS discard list
		for sid in cdna_set - cds_set:
			print(sid, 1, sep="\t", flush=True, file=fout)
