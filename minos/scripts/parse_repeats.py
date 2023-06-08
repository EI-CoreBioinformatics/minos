#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to parse input repeat GFF

Accepted input GFF formats:
- match->match_part
- transcript->exon
"""
import argparse
import csv


def parse_repeats(_in, _out, runid):
    source = runid
    rep_counter = 0
    with open(_out, "w") as out_exons_unstranded:
        for row in csv.reader(open(_in), delimiter="\t"):
            if row and not row[0].startswith("#"):
                if row[2].strip().lower() in {"match_part", "exon"}:
                    rep_counter += 1
                    attrib = dict(
                        item.split("=") for item in row[8].strip(" ;").split(";")
                    )
                    if any(
                        map(
                            lambda x: x is None,
                            (attrib.get("ID"), attrib.get("Parent"),),
                        )
                    ):
                        raise ValueError(
                            "Error: Cannot parse all variables (ID, Parent). Please check entry:\n{}\n".format(
                                "\t".join(row)
                            )
                        )
                    row[1] = source
                    row[2] = "exon"
                    row[6] = "."
                    attrib = "ID={parent}.exon{counter};Parent={parent}".format(
                        parent=attrib["Parent"], counter=rep_counter
                    )
                    print(
                        *row[:8],
                        attrib,
                        sep="\t",
                        file=out_exons_unstranded,
                        flush=True
                    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_gff", type=str)
    ap.add_argument("output_no_strand_exon_gff", type=str)
    ap.add_argument("sample_name", type=str)
    args = ap.parse_args()

    parse_repeats(
        args.input_gff, args.output_no_strand_exon_gff, args.sample_name,
    )


if __name__ == "__main__":
    main()
