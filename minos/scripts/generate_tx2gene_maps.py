import sys
import csv
import re

class Gxf8Parser:
	GFF_REGEX = re.compile("([^;]+)\s*=\s*([^;]+);?")
	GTF_REGEX = re.compile("([^\s;]+)\s+['\"]?([^'\"]+)['\"]?\s*;?")
	GFF_TARGETS = "ID", "Parent"
	GTF_TARGETS = "transcript_id", "gene_id"

	@staticmethod
	def scan_filetype(f):
		for line in open(f):
			if not line.startswith("#"):
				if Gxf8Parser.GFF_REGEX.match(line):
					return Gxf8Parser.GFF_REGEX, Gxf8Parser.GFF_TARGETS
				if Gxf8Parser.GTF_REGEX.match(line):
					return Gxf8Parser.GTF_REGEX, Gxf8Parser.GTF_TARGETS
		return None, None

	def __init__(self, f):
		self.regex, targets = Gxf8Parser.scan_filetype(f)
		if self.regex is None:
			raise ValueError("Cannot determine GTF/GFF type for file {}".format(f))
		self.tid, self.gid = targets

	def parse(self, s):
		return dict((item.group(1).strip(), item.group(2).strip()) for item in self.regex.finditer(s))

	def parse_file(self, _in, runid):
		for row in csv.reader(open(_in), delimiter="\t"):
			if row and not row[0].startswith("#"):
				if row[2].lower() in {"mrna", "ncrna", "transcript"}:
					data = self.parse(row[8])
					yield "{}_{}".format(runid, data[self.tid]), data[self.gid]
