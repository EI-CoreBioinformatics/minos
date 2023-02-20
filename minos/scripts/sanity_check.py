import csv
import os
import argparse
import sys
import csv


class GffValidator:
    gff_header = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]
    valid_feature_types = {
        "gene",
        "mRNA",
        "exon",
        "CDS",
        "five_prime_UTR",
        "three_prime_UTR",
    }

    @staticmethod
    def get_attrib(key, attributes, feature):
        try:
            return attributes[key]
        except:
            raise ValueError(
                "Fatal Error: Cannot parse {} {} field.".format(feature, key)
            )

    def __init__(self, gff_input):
        self.gene_info = dict()
        self.gff_input = gff_input

        for row in csv.DictReader(
            open(gff_input), fieldnames=self.gff_header, delimiter="\t"
        ):
            if not row["seqid"].startswith("#"):
                if row["type"] in {"mRNA", "ncRNA"}:
                    start, end = int(row["start"]), int(row["end"])
                    start, end = (start, end) if start < end else (end, start)
                    attrib = dict(
                        item.split("=")
                        for item in row["attributes"].strip(";").split(";")
                    )

                    gid = GffValidator.get_attrib("Parent", attrib, row["type"]).strip()
                    ginfo = self.gene_info.get(gid, None)
                    if ginfo is None:
                        self.gene_info[gid] = {"start": start, "end": end}
                    else:
                        ginfo["start"] = min(start, ginfo["start"])
                        ginfo["end"] = max(end, ginfo["end"])

    def process(self):
        print("##gff-version 3")
        print_guard = False
        for row in csv.DictReader(
            open(self.gff_input), fieldnames=self.gff_header, delimiter="\t"
        ):
            if not row["seqid"].startswith("#"):
                start, end = int(row["start"]), int(row["end"])
                start, end = (start, end) if start < end else (end, start)
                attrib = dict(
                    item.split("=") for item in row["attributes"].strip(";").split(";")
                )

                if row["type"] == "gene":
                    gid = GffValidator.get_attrib("ID", attrib, row["type"]).strip()
                    ginfo = self.gene_info.get(gid, None)
                    if ginfo is None:
                        raise ValueError(
                            "Error: Something is really wrong. Gene '{}' was not processed before.".format(
                                gid
                            )
                        )

                    start_changed, end_changed = (
                        start != ginfo["start"],
                        end != ginfo["end"],
                    )
                    if start_changed or end_changed:
                        print(
                            "WARN:",
                            gid,
                            "start_changed({} => {})".format(start, ginfo["start"]),
                            "end_changed({} => {})".format(end, ginfo["end"]),
                            sep="\t",
                            file=sys.stderr,
                        )
                    if print_guard:
                        print("###")
                    print_guard = True
                    row["start"], row["end"] = start, end

                print(*row.values(), sep="\t")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gff_input", type=str)

    args = ap.parse_args()

    GffValidator(args.gff_input).process()


if __name__ == "__main__":
    main()

