#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to split and check fasta for BUSCO analysis
"""

import os
from minos import (
    DEFAULT_FAKE_PROT,
    DEFAULT_FAKE_NUC,
)
def split_fasta(fastafile, split_files):
    out = None
    for line in open(fastafile):
        if line.startswith(">"):
            matches = list(
                {runid for runid in split_files if line[1:].startswith(runid + "_")}
            )
            if not matches:
                print("No matching output file for sequence " + line[1:].strip())
                out = None
                continue
            if len(matches) > 1:
                matches = sorted(matches, key=lambda x: len(x), reverse=True)
            out = split_files[matches[0]]
        if out is not None:
            print(line, end="", file=out, flush=True)

    for f in split_files.values():
        f.close()

def copy_file(_in, label, _out):
    with open(_in, 'r') as fh:
        fhdata = fh.read()

    fhdata = fhdata.replace('>', f'>{label}_')

    with open(_out, 'w') as fh:
        fh.write(fhdata)


def check_split_fasta(format, fasta_files):
    for label, path in fasta_files.items():
        fasta = path.name
        if os.path.exists(fasta) and os.stat(fasta).st_size == 0:
            if format == "prot":
                print(f"INFO: Copy fake protein '{DEFAULT_FAKE_PROT}' to empty protein fasta file '{fasta}' to prevent BUSCO analysis failure")
                copy_file(DEFAULT_FAKE_PROT, label, fasta)
            if format == "nuc":
                print(f"INFO: Copy fake nucleotide '{DEFAULT_FAKE_NUC}' to empty nucleotide fasta file '{fasta}' to prevent BUSCO analysis failure")
                copy_file(DEFAULT_FAKE_NUC, label, fasta)