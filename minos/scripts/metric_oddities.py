import sys
import csv
from collections import Counter

# for testing
METRIC_ODDITIES = ['{five_utr_length} >= 10000', '{five_utr_num} >= 5', '{three_utr_length} >= 10000', '{three_utr_num} >= 4', 'not {is_complete}', 'not {has_start_codon}', 'not {has_stop_codon}', '{max_exon_length} >= 10000', '{max_intron_length} >= 500000', '{min_exon_length} <= 5', '{min_intron_length} <= 5', '{selected_cds_fraction} <= 0.3', '{canonical_intron_proportion} != 1', '{non_verified_introns_num} >= 1', 'not {only_non_canonical_splicing}', '{proportion_verified_introns} <= 0.5', '{suspicious_splicing}']


class MetricOddityParser:

	def parse_metrics(self, metric_file, transcript_filter=None):
		for row in csv.DictReader(open(metric_file), delimiter="\t"):
			if transcript_filter is None or row["tid"] in transcript_filter:
				counts = [oddity for oddity in self.oddities if eval(oddity.format(**row))]
				source = row["original_source"]
				self.data.setdefault(
					source,
					Counter({oddity: 0 for oddity in self.oddities})
				).update(counts)
				
	def __init__(self, metric_file, oddities, transcript_filter=None):
		self.data = dict()
		self.oddities = oddities
		self.parse_metrics(metric_file, transcript_filter=transcript_filter)

	def write_table(self, collapse=True, stream=sys.stdout):
		# collapse: True for loci, False for subloci and monoloci
		header = list(self.data.keys()) if not collapse else ["counts"]
		print("metric", *header, sep="\t", file=stream, flush=True)
		for oddity in self.oddities:
			total = 0
			row = list()
			for source, counts in self.data.items():
				total += counts[oddity]
				if not collapse:
					row.append(counts[oddity])
			if collapse:
				row.append(total)
			print(oddity.replace("{", "").replace("}", ""), *row, sep="\t", file=stream, flush=True)
