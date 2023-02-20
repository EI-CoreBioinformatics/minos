import sys
import csv
from collections import Counter


class CollapsedMetricsSummariser:
    def parse_metrics(self, metrics_file):
        self.data = {
            "protein_score": list(),
            "transcript_score": list(),
            "hom_acov_score": list(),
            "protein_score_gene": dict(),
            "transcript_score_gene": dict(),
            "hom_acov_score_gene": dict(),
        }
        for row in csv.DictReader(open(metrics_file), delimiter="\t"):
            for metric in self.data:
                if metric.endswith("_gene"):
                    existing = self.data[metric].get(row["gene"], None)
                    value = float(row[metric])
                    if existing is None:
                        self.data[metric][row["gene"]] = value
                    elif existing != value:
                        raise ValueError(
                            "Ambiguous value found for {metric}: {gene} ({existing}/{value})".format(
                                metric=metric,
                                gene=row["gene"],
                                existing=existing,
                                value=value,
                            )
                        )
                else:
                    self.data[metric].append(float(row[metric]))
        for metric, counts in self.data.items():
            if metric.endswith("_gene"):
                self.data[metric] = list(counts.values())

    def write_summary(self, stream=sys.stdout):
        """
		Average (mean), Count =1, Count >=0.9, Count >=0.75, Count =0, Percentage =1, Percentage >=0.9, Percentage >=0.75, Percentage =0
		"""

        def average(L):
            return sum(L) / len(L)

        limits = ["{} == 1", "{} >= 0.9", "{} >= 0.75", "{} == 0"]
        header = [
            "metric",
            "mean",
            "count = 1",
            "count >= 0.9",
            "count >= 0.75",
            "count = 0",
            "frac = 1",
            "frac >= 0.9",
            "frac >= 0.75",
            "frac = 0",
        ]
        print(*header, sep="\t", file=stream, flush=True)
        for metric, counts in self.data.items():
            row = [metric, average(counts)]
            mlimits = dict()
            for limit in limits:
                mlimits[limit] = sum(eval(limit.format(value)) for value in counts)
            row.extend(mlimits.values())
            row.extend(v / len(counts) for v in mlimits.values())
            print(*row, sep="\t", file=stream, flush=True)

    def __init__(self, metrics_file, stats=None):
        self.parse_metrics(metrics_file)

