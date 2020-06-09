import csv
from collections import Counter

# for testing
METRIC_ODDITIES = ['{five_utr_length} >= 10000', '{five_utr_num} >= 5', '{three_utr_length} >= 10000', '{three_utr_num} >= 4', 'not {is_complete}', 'not {has_start_codon}', 'not {has_stop_codon}', '{max_exon_length} >= 10000', '{max_intron_length} >= 500000', '{min_exon_length} <= 5', '{min_intron_length} <= 5', '{selected_cds_fraction} <= 0.3', '{canonical_intron_proportion} != 1', '{non_verified_introns_num} >= 1', 'not {only_non_canonical_splicing}', '{proportion_verified_introns} <= 0.5', '{suspicious_splicing}']


class MetricOddityParser:
	def __init__(self, metric_file, oddities, gene_filter=None):
		self.table = Counter({oddity: 0 for oddity in oddities})
		self.metric_file = metric_file
		self.gene_filter = gene_filter
	def run(self):
		for row in csv.DictReader(open(self.metric_file), delimiter="\t"):
			if gene_filter is None or row["tid"] in gene_filter:
				counts = [oddity for oddity in self.table if eval(oddity.format(**row))]
				self.table.update(counts)

		for oddity, count in self.table.items():
			print(oddity.replace("{", "").replace("}", ""), count, sep="\t")


				
