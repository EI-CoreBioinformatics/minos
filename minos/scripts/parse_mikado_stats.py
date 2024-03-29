import csv
import argparse
import os
import sys

STATS_ORDER = [
    "Number of genes",
    "Number of Transcripts",
    "Transcripts per gene",
    "Number of monoexonic genes",
    "Monoexonic transcripts",
    "Transcript mean size cDNA (bp)",
    "Transcript median size cDNA (bp)",
    "Min cDNA",
    "Max cDNA",
    "Total exons",
    "Exons per transcript",
    "Exon mean size (bp)",
    "CDS mean size (bp)",
    "Transcript mean size CDS (bp)",
    "Transcript median size CDS (bp)",
    "Min CDS",
    "Max CDS",
    "Intron mean size (bp)",
    "5'UTR mean size (bp)",
    "3'UTR mean size (bp)",
]


def format_number(x):
    return "{:.02f}".format(float(x.replace(",", "")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stats_tsv", type=str)
    args = ap.parse_args()

    stats_info = dict()
    # Stat    Total   Average Mode    Min     1%      5%      10%     25%     Median  75%     90%     95%     99%     Max
    for row in csv.DictReader(open(args.stats_tsv), delimiter="\t"):
        try:
            stat = row["Stat"].strip()
        except:
            continue

        if stat.lower() in {
            "number of genes",
            "transcripts per gene",
            "number of monoexonic genes",
            "monoexonic transcripts",
            "exons per transcript",
        }:
            stats_info[stat] = format_number(row["Total"])
            if stat.lower() == "transcripts per gene":
                stats_info["Number of Transcripts"] = stats_info[stat]
                stats_info[stat] = format_number(row["Average"])
            elif stat.lower() == "exons per transcript":
                stats_info["Total exons"] = stats_info[stat]
                stats_info[stat] = format_number(row["Average"])

        elif stat.lower() in {
            "exon lengths",
            "cds exon lengths",
            "intron lengths",
            "5'utr length",
            "3'utr length",
        }:
            stats_info["{} mean size (bp)".format(stat.split(" ")[0])] = format_number(
                row["Average"]
            )

        elif stat.lower() in {"cdna lengths", "cds lengths"}:
            stats_info[
                "Transcript mean size {} (bp)".format(stat.split(" ")[0])
            ] = format_number(row["Average"])
            stats_info[
                "Transcript median size {} (bp)".format(stat.split(" ")[0])
            ] = format_number(row["Median"])
            stats_info["Min {}".format(stat.split(" ")[0])] = format_number(row["Min"])
            stats_info["Max {}".format(stat.split(" ")[0])] = format_number(row["Max"])

    for stat in STATS_ORDER:
        print(stat, stats_info[stat.replace("cDNA", "CDNA")], sep="\t")


if __name__ == "__main__":
    main()
