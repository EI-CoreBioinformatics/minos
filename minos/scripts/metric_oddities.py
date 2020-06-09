import sys
import csv
from collections import Counter

# for testing
METRIC_ODDITIES = ['{five_utr_length} >= 10000', '{five_utr_num} >= 5', '{three_utr_length} >= 10000', '{three_utr_num} >= 4', 'not {is_complete}', 'not {has_start_codon}', 'not {has_stop_codon}', '{max_exon_length} >= 10000', '{max_intron_length} >= 500000', '{min_exon_length} <= 5', '{min_intron_length} <= 5', '{selected_cds_fraction} <= 0.3', '{canonical_intron_proportion} != 1', '{non_verified_introns_num} >= 1', 'not {only_non_canonical_splicing}', '{proportion_verified_introns} <= 0.5', '{suspicious_splicing}']


class MetricOddityParser:

	def parse_metrics(self, metric_file, transcript_filter=None):
		data = dict()
		for row in csv.DictReader(open(metric_file), delimiter="\t"):
			if transcript_filter is None or row["tid"] in transcript_filter:
				counts = [oddity for oddity in self.oddities if eval(oddity.format(**row))]
				source = row["original_source"]
				data.setdefault(
					source,
					Counter({oddity: 0 for oddity in self.oddities})
				).update(counts)
		return data

	def calc_stats(self, data, collapse=True):
		stats = dict()
		for oddity in self.oddities:
			total = 0
			row = list()
			for source, counts in data.items():
				total += counts[oddity]
				if not collapse:
					row.append(counts[oddity])
			if collapse:
				row.append(total)
			stats[oddity.replace("{", "").replace("}", "")] = row
		return stats

	def __init__(self, final_metric_file, subloci_metric_file, monoloci_metric_file, oddities, transcript_filter=None):
		self.oddities = oddities
		self.data = self.calc_stats(self.parse_metrics(final_metric_file, transcript_filter=transcript_filter))
		sub_data = self.parse_metrics(subloci_metric_file)
		mono_data = self.parse_metrics(monoloci_metric_file)
		combined = dict()
		for tset in [sub_data, mono_data]:
			for source, counts in tset.items():
				combined.setdefault(source, Counter()).update(counts)
		for metric, values in self.calc_stats(combined, collapse=False).items():
			self.data[metric].extend(values)
		self.header = ["metric", "final_set"] + list(combined.keys())

	def write_table(self, stream=sys.stdout):
		print(*self.header, sep="\t", file=stream, flush=True)
		for metric, values in self.data.items():
			print(metric, *values, sep="\t", file=stream, flush=True)
