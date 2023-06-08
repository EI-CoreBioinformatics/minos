#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to collapse metrics
"""

import argparse
import sys
import csv
import os
import yaml
from eicore.snakemake_helper import DEFAULT_CONFIG_FILE

run_config = dict()
minos_config = str()
# use minos config from ENV by default or use default minos config
if "DEFAULT_MINOS_CONFIG_FILE" in os.environ:
    minos_config = os.path.abspath(os.environ["DEFAULT_MINOS_CONFIG_FILE"])
    if os.path.exists(minos_config):
        run_config = yaml.load(open(minos_config), Loader=yaml.SafeLoader)
    else:
        print(
            f"WARNING: Could not find file associated with ENV DEFAULT_MINOS_CONFIG_FILE: '{minos_config}'. Going to use the default minos config instead {DEFAULT_CONFIG_FILE}",
            file=sys.stderr,
        )
        minos_config = DEFAULT_CONFIG_FILE
        run_config = yaml.load(open(DEFAULT_CONFIG_FILE), Loader=yaml.SafeLoader)

else:
    minos_config = DEFAULT_CONFIG_FILE
    run_config = yaml.load(open(DEFAULT_CONFIG_FILE), Loader=yaml.SafeLoader)

print(f"INFO: Minos config used: '{minos_config}'", file=sys.stderr)

SCORES = [
    "protein_score",
    "transcript_score",
    "hom_qcov_score",
    "hom_tcov_score",
    "hom_acov_score",
    "te_score",
    "cpc_score",
    "expression_score",
    "classification",
    "busco_score",
]

#  3 + 9 * 2 + 5 = 26
HEADER = (
    ["#transcript", "gene", "alias"]
    + SCORES
    + [s + "_gene" for s in SCORES]
    + ["confidence", "repeat_associated", "biotype", "discard", "region"]
)


def check_expression(expression, values):
    def cmp_score(a, b, op):
        if op == "eq":
            return a == b
        if op == "lt":
            return a < b
        if op == "gt":
            return a > b
        if op == "ge":
            return a >= b
        raise ValueError("Invalid check: {} {} {}".format(a, b, op))

    if type(expression[0]) is str:
        a, b, op = expression
        expression = (values[a], b, op)
        return cmp_score(*expression)
    operator = expression[0]
    return operator(check_expression(exp, values) for exp in expression[1:])


class TranscriptData(dict):
    def __init__(self, gff):
        try:
            for row in csv.reader(open(gff), delimiter="\t"):
                if not row[0].startswith("#"):
                    ftype = row[2].strip().lower()
                    if ftype in {"mrna", "ncrna"}:
                        # TODO: The attribute strip to be made more comprehensive
                        attrib = dict(
                            item.split("=") for item in row[8].strip("; ").split(";")
                        )

                        if any(
                            map(
                                lambda x: x is None,
                                (
                                    attrib.get("ID"),
                                    attrib.get("Parent"),
                                    attrib.get("Name"),
                                ),
                            )
                        ):
                            raise ValueError(
                                "Error: Cannot parse all variables (ID, Parent, Name). Please check entry:\n{}\n".format(
                                    "\t".join(row)
                                )
                            )

                        if self.get(attrib["Name"]) is not None:
                            raise ValueError(
                                "Error: Potential duplicate entry. Transcript '{}' already processed. Please check.\n{}\n".format(
                                    attrib["Name"], "\t".join(row)
                                )
                            )

                        start, end = map(int, row[3:5])
                        start, end = (start, end) if start < end else (end, start)

                        self[attrib["Name"]] = {
                            "id": attrib["ID"],
                            "parent": attrib["Parent"],
                            "alias": attrib["Name"],
                            "region": "{}:{}..{}".format(row[0], start, end),
                            "type": ftype,
                        }
        except FileNotFoundError:
            raise FileNotFoundError("Cannot find input gff at " + gff)


class ExpressionData(dict):
    def read_kallisto(self, tsv):
        for row in csv.DictReader(open(tsv), delimiter="\t"):
            prev_tpm = self.get(row["target_id"], 0.0)
            try:
                tpm = float(row["tpm"])
            except:
                tpm = 0.0
            self[row["target_id"]] = max(prev_tpm, tpm)

    def __init__(self, *expression_data_files):
        files_imported = 0
        for f in expression_data_files:
            try:
                self.read_kallisto(f)
            except FileNotFoundError:
                print("Warning: could not find kallisto data file at " + f)
            else:
                files_imported += 1

        if files_imported == 0:
            print(
                "Warning: No kallisto data processed. Assigning kallisto scores to 0.0"
            )


class TranscriptScores:
    def __init__(self, metrics_info, expression_score, short_cds, **metrics):
        self.tid = metrics["tid"]
        self.classification = 0  ## cschu 20200203: issue9: disabled until full-lengther replacement implemented
        self.short_cds = short_cds
        self.expression_score = expression_score

        self.protein_score = 0
        self.protein_score = max(
            [0]
            + [
                float(metrics[m + "_aF1"])
                for m in metrics_info.get("mikado.protein", set())
            ]
        )
        self.transcript_score = max(
            [0]
            + [
                float(metrics[m + "_aF1"])
                for m in metrics_info.get("mikado.transcript", set())
            ]
        )
        self.hom_qcov_score = max(
            [0]
            + [float(metrics[m + "_qCov"]) for m in metrics_info.get("blast", set())]
        )
        self.hom_tcov_score = max(
            [0]
            + [float(metrics[m + "_tCov"]) for m in metrics_info.get("blast", set())]
        )
        self.hom_acov_score = (self.hom_qcov_score + self.hom_tcov_score) / 2.0
        self.te_score = max(
            [0]
            + [float(metrics[m + "_cov"]) for m in metrics_info.get("repeat", set())]
        )
        # 0.0  # !TODO, # we get the highest for the te and when we compute for the gene we take the lowest downstream
        self.cpc_score = float(metrics["cpc"])
        self.busco_score = float(metrics["busco_proteins"])


class MetricCollapser:
    def read_metrics_info(self, tsv):
        self.metrics_info = dict()
        try:
            for row in csv.reader(open(tsv), delimiter="\t"):
                metric_type, metric_id = row[0:2]
                self.metrics_info.setdefault(metric_type, set()).add(metric_id)
        except FileNotFoundError:
            raise FileNotFoundError("Cannot find metrics info at " + tsv)

    def read_metrics(self, matrix):
        try:
            for metrics in csv.DictReader(open(matrix), delimiter="\t"):
                tdata = self.transcripts_data.get(metrics["tid"])
                # print(metrics["tid"], metrics["tid"] in self.transcripts_data)
                if tdata is not None:
                    mdata = self.model_info.get(metrics["tid"])
                    if mdata is not None:
                        raise ValueError(
                            "Error: Potential duplicate entry. Transcript '{}' already processed. Please check.\n{}\n".format(
                                metrics["tid"], "\t".join(metrics)
                            )
                        )

                    kallisto_score = self.expression_data.get(tdata["id"], 0.0)
                    if kallisto_score is None:
                        raise ValueError(
                            "Error: Could not extract tpm data for transcript {} ({})".format(
                                metrics["tid"], tdata["id"]
                            )
                        )

                    self.model_info[metrics["tid"]] = dict(tdata)
                    gid = self.model_info[metrics["tid"]]["gene"] = tdata["parent"]
                    del self.model_info[metrics["tid"]]["parent"]
                    tmetrics = TranscriptScores(
                        self.metrics_info,
                        kallisto_score,
                        self.short_cds[tdata["id"]],
                        **metrics,
                    ).__dict__
                    self.model_info[metrics["tid"]].update(tmetrics)

                    #  get the highest metrics value for gene
                    # except for te_score
                    ginfo = self.gene_info.get(gid, dict())
                    if not ginfo:
                        self.gene_info[gid] = dict(tmetrics)
                    else:
                        for k, v in ginfo.items():
                            cmp_f = max if k not in {"te_score", "short_cds"} else min
                            self.gene_info[gid][k] = cmp_f(v, tmetrics[k])

        except FileNotFoundError:
            raise FileNotFoundError("Cannot find metrics matrix at " + matrix)

    def __init__(self, gff, metrics_info, metrics_matrix, short_cds, expression_data):
        self.read_metrics_info(metrics_info)
        self.transcripts_data = TranscriptData(gff)
        self.expression_data = ExpressionData(*expression_data)
        self.short_cds = dict(
            (row[0], int(row[1])) for row in csv.reader(open(short_cds), delimiter="\t")
        )

        self.model_info, self.gene_info = dict(), dict()
        self.read_metrics(metrics_matrix)

    def write_scores(self, checks, stream=sys.stdout):
        print(*HEADER, sep="\t", file=stream)

        for tid, tinfo in sorted(self.model_info.items(), key=lambda x: x[0]):
            row = [tinfo["id"], tinfo["gene"], tinfo["alias"]]
            row.extend(tinfo.get(score, ".") for score in SCORES)
            gid = tinfo["gene"]
            gene_scores = list(
                self.gene_info.get(gid, dict()).get(score, ".") for score in SCORES
            )
            if any(score == "." for score in gene_scores):
                raise ValueError(
                    "Error: Cannot find all gene scores:\n{}".format(
                        "\n".join(zip(SCORES, gene_scores))
                    )
                )
            row.extend(gene_scores)

            # classify protein coding and ncRNA models

            # classify ncRNA models
            if tinfo["type"].lower() == "ncrna":
                discard = False
                ncrna_gene = eval(checks["ncrna_gene"].format(**self.gene_info[gid]))
                biotype = "ncrna_gene"
                # discard any ncRNA not meeting the collapse_metrics_thresholds
                if not ncrna_gene:
                    discard = True
                repeat_associated = eval(
                    checks["repeat_associated"].format(**self.gene_info[gid])
                )
                if repeat_associated:
                    biotype = "transposable_element_gene"
                # standard for ncRNA for now
                high_confidence = False
                row.extend(
                    [
                        "High" if high_confidence else "Low",
                        str(repeat_associated),
                        biotype,
                        str(discard),
                        tinfo["region"],
                    ]
                )
            # classify protein coding models
            else:
                biotype = "protein_coding_gene"
                repeat_associated = eval(
                    checks["repeat_associated"].format(**self.gene_info[gid])
                )
                if repeat_associated:
                    biotype = "transposable_element_gene"
                elif eval(checks["predicted_gene"].format(**self.gene_info[gid])):
                    biotype = "predicted_gene"

                high_confidence = eval(
                    checks["hi_confidence"].format(**self.gene_info[gid])
                )
                discard = eval(checks["discard_coding"].format(**self.gene_info[gid]))

                row.extend(
                    [
                        "High" if high_confidence else "Low",
                        str(repeat_associated),
                        biotype,
                        str(discard),
                        tinfo["region"],
                    ]
                )

            print(*row, sep="\t", file=stream)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_gff", type=str)
    ap.add_argument("metrics_matrix", type=str)
    ap.add_argument("metrics_info", type=str)
    ap.add_argument("short_cds", type=str)
    ap.add_argument("kallisto_tpm", nargs="*")
    args = ap.parse_args()

    mc = MetricCollapser(
        args.input_gff,
        args.metrics_info,
        args.metrics_matrix,
        args.short_cds,
        args.kallisto_tpm,
    )
    mc.write_scores(run_config["collapse_metrics_thresholds"])


if __name__ == "__main__":
    main()

