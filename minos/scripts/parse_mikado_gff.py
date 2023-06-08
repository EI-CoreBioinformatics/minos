import sys
import argparse
import csv
import collections


def get_alias(tid, attrib):

    err_msg = "# {} attribute not defined for transcript {}. {}"
    resolution = "Checking {} instead ... "
    success_msg = "# Got {}, using it as Name attribute."

    try:
        alias = attrib["alias"]
    except:
        print(
            err_msg.format("alias", tid, resolution.format("target")),
            end="",
            file=sys.stderr,
        )
        try:
            alias = attrib["target"]
            print(success_msg.format("target"), file=sys.stderr)
        except:
            print(
                err_msg.format("target", tid, resolution.format("prev_parent")),
                end="",
                file=sys.stderr,
            )
            try:
                alias = attrib["prev_parent"]
                print(success_msg.format("prev_parent"), file=sys.stderr)
            except:
                print(
                    err_msg.format("prev_parent", tid, "Falling back to ID ... "),
                    end="",
                    file=sys.stderr,
                )
                alias = tid

    return alias[0] if type(alias) is list else alias


class GffAttributes(dict):
    def __init__(self, attribs):
        for item in attribs.split(";"):
            try:
                k, v = item.split("=")
            except:
                pass
            self.setdefault(k, list()).append(v)

    def __str__(self):
        return ";".join("{}={}".format(k, v[0]) for k, v in self.items())


class GffRecord:
    __fields = [
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

    def __init__(self, data):
        try:
            for k, v in zip(self.__fields[:-1], data[:-1]):
                setattr(self, k, v)
        except:
            raise ValueError(
                "Error when parsing gff line:\n{}\n".format("\t".join(data))
            )

        self.attributes = GffAttributes(data[-1])

    def is_feature(self, feature):
        return self.type.lower() == feature.lower()

    def get_attribute(self, attr_id):
        try:
            return self.attributes[attr_id]
        except KeyError:
            raise KeyError(
                "Error: Could not find attribute {} in\n{}\n".format(attr_id, str(self))
            )

    def normalize_feature(self):
        norm_map = {
            "sublocus": "gene",
            "ncrna_gene": "ncRNA_gene",
            "mrna": "mRNA",
            "ncrna": "ncRNA",
            "transcript": "mRNA",
        }

        if self.type.lower() == "sublocus":
            self.attributes["superlocus"] = self.get("Parent", ".")
            try:
                del self["Parent"]
            except:
                pass

        try:
            self.type = norm_map[self.type.lower()]
        except:
            raise ValueError(
                "Error: Invalid feature {} type in \n{}\n".format(self.type, str(self))
            )

    def __str__(self):
        sentinel = "###\n" if self.is_feature("gene") else ""
        return sentinel + "\t".join(str(getattr(self, attr)) for attr in self.__fields)


class MikadoLociParser:
    @staticmethod
    def parse(gff_file, source="."):
        for row in csv.reader(open(gff_file), delimiter="\t"):
            if row[0]:
                if row[0].startswith("#"):
                    print(*row, sep="\t")
                else:
                    gff_rec = GffRecord(row)
                    gff_rec.source = source

                    if gff_rec.is_feature("superlocus"):
                        continue
                    elif gff_rec.is_feature("sublocus"):
                        gff_rec.normalize_feature()
                    elif gff_rec.is_feature("ncRNA_gene"):
                        gff_rec.normalize_feature()
                    elif any(
                        gff_rec.is_feature(feat)
                        for feat in ("mRNA", "ncRNA", "transcript")
                    ):
                        gff_rec.normalize_feature()
                        tid = gff_rec.get_attribute("ID")[0]

                        preserve_attribs = ("ID", "Parent", "Name", "alias", "Note")
                        new_attrib = {
                            "ID": tid,
                            "Parent": gff_rec.get_attribute("Parent")[0],
                            "Name": get_alias(tid, gff_rec.attributes),
                            "Note": [tid],
                        }
                        new_attrib["Note"].extend(
                            "{}:{}".format(k, v[0])
                            for k, v in gff_rec.attributes.items()
                            if k not in preserve_attribs and k.lower() != "note"
                        )
                        new_attrib["Note"].extend(
                            gff_rec.attributes.get("note", [""])[0].split("|")[1:]
                        )
                        new_attrib["Note"] = ",".join(new_attrib["Note"])
                        gff_rec.attributes = ";".join(
                            "{}={}".format(k, v) for k, v in new_attrib.items()
                        )

                    print(gff_rec)


def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("gff", type=str)
    ap.add_argument("--source", type=str, default="Mikado_annotation_run2")
    args = ap.parse_args()

    MikadoLociParser.parse(args.gff, source=args.source)


if __name__ == "__main__":
    main()

