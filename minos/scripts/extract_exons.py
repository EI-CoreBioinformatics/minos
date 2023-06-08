#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to extract exon from GTF file and write output just exons GFF format
"""
import argparse
import csv
import re


def extract_exons(_in, _out):
    with open(_out, "w") as exons_out:
        exon = 1
        for row in csv.reader(open(_in), delimiter="\t"):
            if row and not row[0].startswith("#") and row[2].lower() == "exon":
                tid = str()
                try:
                    tid = re.search('transcript_id "(.+?)"', row[8]).group(1)
                except AttributeError as err:
                    raise AttributeError(
                        "Error: {err}.\nCannot parse all variables (transcript_id). Please check entry:\n{line}\n".format(
                            err=err, line="\t".join(row)
                        )
                    )
                row[8] = f"ID={tid}.exon{exon};Parent={tid}"
                exon += 1
                print(*row, sep="\t", flush=True, file=exons_out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_gtf", type=str)
    ap.add_argument("output_exon_gff", type=str)
    args = ap.parse_args()

    extract_exons(args.input_gtf, args.output_exon_gff)


if __name__ == "__main__":
    main()
