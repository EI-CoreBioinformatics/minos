#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to clean fasta header
"""

# import libraries
import sys
import argparse


class CleanFastaHeader(object):
    def __init__(self, args):
        self.args = args

    def clean_fasta_header(self):
        for line in self.args.fasta:
            line = line.rstrip()
            if line.startswith(">"):
                x = line.split(" ")
                header = x[0]
                if self.args.add_fields:
                    cols = (y.strip()
                            for y in self.args.add_fields.split(","))
                    header = header + " " + \
                        " ".join(
                            i for j in cols for i in x if i.lower().startswith(j.lower()))
                if self.args.add_columns:
                    cols = (int(y.strip()) -
                            1 for y in self.args.add_columns.split(","))
                    header = header + " " + " ".join(x[z] for z in cols)
                print(f"{header}")
            else:
                print(line)

    def run(self):
        self.clean_fasta_header()


def main():
    parser = argparse.ArgumentParser(
        description="Script to clean fasta header")
    parser.add_argument("fasta", nargs='?', type=argparse.FileType(
        'r'), default=sys.stdin, help="Provide fasta file (as a file or stdin)")
    parser.add_argument("--add_fields", nargs='?',
                        help="Add these additional fields to header (comma delimited). For example, giving 'Name,Note' (with quotes) will add those fields if fasta header has fields that starts with them. The fasta headers are split by space. Please avoid using '>' as one of the field. Case-insensitive.")
    parser.add_argument("--add_columns", nargs='?',
                        help="Add these additional columns to header (comma delimited). For example, giving '2,3,4' (with quotes) will add those fields to the fasta header. The fasta headers are split by space. Column 1 is first column, which is the fasta header, so please avoid using that")
    args = parser.parse_args()

    CleanFastaHeader(args).run()


if __name__ == "__main__":
    main()
