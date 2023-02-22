import argparse
import sys
import csv


class TranscriptDataValidator:
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
        "ncRNA_gene",
        "mRNA",
        "ncRNA",
        "exon",
        "CDS",
        "five_prime_UTR",
        "three_prime_UTR",
    }

    def __init__(self, gff):
        self.genes_info = dict()
        self.transcripts_info = dict()
        self.transcript_gene_info = dict()
        self.exons_info = dict()
        self.exon_coords = dict()
        self.cds_info = dict()
        self.cds_coords = dict()
        self.cds_length = dict()
        self.utr5_info = dict()
        self.utr3_info = dict()

        self.__read_gff(gff)

    def __read_gff(self, gff):
        def get_attrib(key, attributes, feature):
            try:
                return attributes[key]
            except:
                raise ValueError(
                    "Fatal Error: Cannot parse {} {} field.".format(feature, key)
                )

        try:
            with open(gff) as gff_in:
                for row in csv.DictReader(
                    open(gff), fieldnames=self.gff_header, delimiter="\t"
                ):
                    if not row["seqid"].startswith("#"):
                        if row["type"] in self.valid_feature_types:
                            start, end = int(row["start"]), int(row["end"])
                            if start > end:
                                print(
                                    "WARN: start coordinate is greater than end coordinate [{} > {}], swapping them"
                                )
                                row["start"], row["end"] = end, start
                            else:
                                row["start"], row["end"] = start, end

                            if row["strand"] not in {"+", "-"}:
                                # if a strand is a dot, then assume it is a positive one
                                print(
                                    "WARN: strand is not defined as + or -, assuming it is +, plus strand"
                                )
                                row["strand"] = "+"

                            # make sure phase is an integer
                            row["phase"] = (
                                int(row["phase"])
                                if row["phase"] in {"0", "1", "2"}
                                else 0
                            )

                            attrib = dict(
                                item.split("=")
                                for item in row["attributes"].strip(";").split(";")
                            )

                            if row["type"] == "gene":
                                gid = get_attrib("ID", attrib, row["type"])

                                if self.genes_info.get(gid, None) is not None:
                                    raise ValueError(
                                        "Fatal error: Duplicate gene ID '{}' in input file {}".format(
                                            gid, gff
                                        )
                                    )

                                self.genes_info[gid] = {
                                    "start": row["start"],
                                    "end": row["end"],
                                }

                            elif row["type"] == "mRNA":
                                tid = get_attrib("ID", attrib, row["type"])
                                gid = get_attrib("Parent", attrib, row["type"])

                                if self.transcripts_info.get(tid, None) is not None:
                                    raise ValueError(
                                        "Fatal error: Duplicate mRNA ID '{}' in input file {}".format(
                                            tid, gff
                                        )
                                    )

                                self.transcripts_info[tid] = {
                                    "start": row["start"],
                                    "end": row["end"],
                                    "strand": row["strand"],
                                    "gene": gid,
                                }

                                ginfo = self.transcript_gene_info.get(gid, None)
                                if ginfo is None:
                                    self.transcript_gene_info[gid] = {
                                        "start": row["start"],
                                        "end": row["end"],
                                    }
                                else:
                                    ginfo.update(
                                        {
                                            "start": min(ginfo["start"], row["start"]),
                                            "end": max(ginfo["end"], row["end"]),
                                        }
                                    )

                            elif row["type"] == "exon":
                                eid = get_attrib("Parent", attrib, row["type"])
                                self.exons_info.setdefault(eid, list()).append(
                                    (row["start"], row["end"])
                                )

                                einfo = self.exon_coords.get(eid, None)
                                if einfo is None:
                                    self.exon_coords[eid] = {
                                        "start": row["start"],
                                        "end": row["end"],
                                    }
                                else:
                                    einfo.update(
                                        {
                                            "start": min(einfo["start"], row["start"]),
                                            "end": max(einfo["end"], row["end"]),
                                        }
                                    )

                            elif row["type"] == "CDS":
                                cid = get_attrib("Parent", attrib, row["type"])
                                self.cds_info.setdefault(cid, list()).append(
                                    (row["start"], row["end"])
                                )

                                cinfo = self.cds_coords.get(cid, None)
                                if cinfo is None:
                                    self.cds_coords[cid] = {
                                        "start": row["start"],
                                        "end": row["end"],
                                    }
                                else:
                                    cinfo.update(
                                        {
                                            "start": min(cinfo["start"], row["start"]),
                                            "end": max(cinfo["end"], row["end"]),
                                        }
                                    )

                                # deduct the phase from the end
                                # and add to the hash cds_length_info when calculating CDS length
                                self.cds_length.setdefault(cid, list()).append(
                                    (row["start"], row["end"] - row["phase"])
                                )

                            elif (
                                row["type"] == "five_prime_UTR"
                                or row["type"] == "three_prime_UTR"
                            ):
                                uid = get_attrib(
                                    "Parent", attrib, row["type"].split("_")[-1]
                                )
                                if row["type"] == "five_prime_UTR":
                                    self.utr5_info.setdefault(uid, list()).append(
                                        (row["start"], row["end"])
                                    )
                                if row["type"] == "three_prime_UTR":
                                    self.utr3_info.setdefault(uid, list()).append(
                                        (row["start"], row["end"])
                                    )

        except FileNotFoundError:
            raise FileNotFoundError("Error: Cannot find input gff at " + gff)

        for d in [
            self.utr5_info,
            self.utr3_info,
            self.exons_info,
            self.cds_info,
            self.cds_length,
        ]:
            for v in d.values():
                v.sort()

    def validate_gene_spans(self):
        for gid, mginfo in sorted(
            self.transcript_gene_info.items(), key=lambda x: x[0]
        ):
            ginfo = self.genes_info.get(gid, None)
            if ginfo is None:
                print("WARN: Gene id '{}' not found.".format(gid))
            if ginfo["start"] != mginfo["start"]:
                print(
                    "WARN: Gene start is not consistent to all mRNA spans for gene '{}' (actual={}, computed={})".format(
                        gid, ginfo["start"], mginfo["start"]
                    )
                )
            if ginfo["end"] != mginfo["end"]:
                print(
                    "WARN: Gene end is not consistent to all mRNA spans for gene '{}' (actual={}, computed={})".format(
                        gid, ginfo["end"], mginfo["end"]
                    )
                )

    def validate_transcripts(self):
        for tid, tinfo in sorted(self.transcripts_info.items(), key=lambda x: x[0]):
            actual_cds_start, actual_cds_end = sorted(
                self.cds_coords.get(tid, {"start": 0, "end": 0}).values()
            )

            for f in ["exon", "cds", "utr5", "utr3"]:
                for v in ["count", "len", "start", "end"] + (
                    ["complete_match", "partial_match"] if f != "exon" else []
                ):
                    tinfo["{}_{}".format(f, v)] = 0
                    if f == "cds" and v == "len":
                        tinfo["cds_len_standard"] = 0

            for start, end in self.exons_info.get(tid, list()):
                span = end - start + 1
                tinfo["exon_count"] += 1
                tinfo["exon_len"] += span
                if tinfo["exon_start"] == 0 or start < tinfo["exon_start"]:
                    tinfo["exon_start"] = start
                if tinfo["exon_end"] == 0 or end > tinfo["exon_end"]:
                    tinfo["exon_end"] = end

            for start, end in self.cds_info.get(tid, list()):
                span = end - start + 1
                tinfo["cds_count"] += 1
                tinfo["cds_len_standard"] += span

                for e_start, e_end in self.exons_info.get(tid, list()):
                    if e_start == start and e_end == end:
                        tinfo["cds_complete_match"] += 1
                    elif (
                        (e_start < start and end == e_end)
                        or (e_start == start and end < e_end)
                        or (e_start < start and end < e_end)
                    ):
                        tinfo["cds_partial_match"] += 1

                total_cds_count = (
                    tinfo["cds_complete_match"] + tinfo["cds_partial_match"]
                )
                if total_cds_count != tinfo["cds_count"]:
                    print(
                        "WARN: Not all CDS exons for transcript '{}' are covered within exon.".format(
                            tid
                        )
                    )

            for start, end in self.cds_length.get(tid, list()):
                span = end - start + 1
                tinfo["cds_len"] += span

            for utr, utr_info in [(5, self.utr5_info), (3, self.utr3_info)]:
                for start, end in utr_info.get(tid, list()):
                    span = end - start + 1
                    tinfo["utr{}_count".format(utr)] += 1
                    tinfo["utr{}_len".format(utr)] += span

                    for e_start, e_end in self.exons_info.get(tid, list()):
                        if e_start == start and end == e_end:
                            tinfo["utr{}_complete_match".format(utr)] += 1
                        elif (
                            utr == 5
                            and tinfo["strand"] == "+"
                            and e_start == start
                            and end < e_end
                        ):
                            tinfo["utr{}_partial_match".format(utr)] += 1
                            if end != actual_cds_start - 1:
                                print(
                                    "WARN: {}'-UTR start position not at expected start for transcript '{}'".format(
                                        utr, tid
                                    )
                                )
                        elif (
                            utr == 5
                            and tinfo["strand"] == "-"
                            and e_start < start
                            and end == e_end
                        ):
                            tinfo["utr{}_partial_match".format(utr)] += 1
                            if start != actual_cds_end + 1:
                                print(
                                    "WARN: {}'-UTR start position not at expected start for transcript '{}'".format(
                                        utr, tid
                                    )
                                )
                        elif (
                            utr == 3
                            and tinfo["strand"] == "+"
                            and e_start < start
                            and end == e_end
                        ):
                            tinfo["utr{}_partial_match".format(utr)] += 1
                            if start != actual_cds_end + 1:
                                print(
                                    "WARN: {}'-UTR start position not at expected start for transcript '{}'".format(
                                        utr, tid
                                    )
                                )
                        elif (
                            utr == 3
                            and tinfo["strand"] == "-"
                            and e_start == start
                            and end < e_end
                        ):
                            tinfo["utr{}_partial_match".format(utr)] += 1
                            if end != actual_cds_start - 1:
                                print(
                                    "WARN: {}'-UTR start position not at expected start for transcript '{}'".format(
                                        utr, tid
                                    )
                                )

                    total_utr_count = (
                        tinfo["utr{}_complete_match".format(utr)]
                        + tinfo["utr{}_partial_match".format(utr)]
                    )
                    if total_utr_count != tinfo["utr{}_count".format(utr)]:
                        print(
                            "WARN: Not all {}'-UTR exons for transcript '{}' are covered within exon.".format(
                                utr, tid
                            )
                        )

            # check to make sure that we have exon features for transcripts with CDS feature.
            if tinfo["exon_len"] == 0 and (
                tinfo["cds_len"] > 0 or tinfo["cds_len_standard"] > 0
            ):
                print(
                    "WARN: Transcript '{}' does not have exon feature, but has CDS feature".format(
                        tid
                    )
                )
                cds_cdna_ratio = 0
            else:
                cds_cdna_ratio = tinfo["cds_len"] / tinfo["exon_len"]
                if tinfo["cds_len_standard"] > tinfo["exon_len"]:
                    print(
                        "WARN: Transcript '{}' CDS length is greater than exon length".format(
                            tid
                        )
                    )

            if tinfo["start"] != tinfo["exon_start"]:
                print(
                    "WARN: Transcript start is not consistent with all exon spans for transcript '{}'.".format(
                        tid
                    )
                )
            if tinfo["end"] != tinfo["exon_end"]:
                print(
                    "WARN: Transcript end is not consistent with all exon spans for transcript '{}'.".format(
                        tid
                    )
                )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_gff", type=str)
    args = ap.parse_args()

    tdv = TranscriptDataValidator(args.input_gff)
    tdv.validate_gene_spans()
    tdv.validate_transcripts()


if __name__ == "__main__":
    main()

