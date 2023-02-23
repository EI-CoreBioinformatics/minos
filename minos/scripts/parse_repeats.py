import csv
import re

# awk '$3 == "match_part"' {input[0]} | sed -e 's/\\tmatch_part\\t/\\texon\\t/' -e 's/\\t[+-]\\t/\\t.\\t/' > {output[0]}
def parse_repeatmasker(_in, _out, _out2, runid):
    source = tag = runid
    name_tag, rep_counter = "RM", 0
    with open(_out, "w") as out_exons, open(_out2, "w") as out_exons_unstranded:
        for row in csv.reader(open(_in), delimiter="\t"):
            if row and not row[0].startswith("#") and row[1] == "RepeatMasker":
                rep_counter += 1
                try:
                    target = re.search('Target\s+"([^"]+)', row[8]).group(1)
                except:
                    try:
                        target = re.search("Target\s*=\s*([^;]+)", row[8]).group(1)
                    except:
                        raise ValueError("Cannot parse Target from " + row[8])
                target = re.sub("\s+", "_", target)
                note = ";Note={}".format(target)
                try:
                    name = re.sub(
                        "\s+", "_", re.search("Name\s*=\s*([^;]+)", row[8]).group(1)
                    )
                except:
                    name, note = target, ""
                row[1] = source
                row[5] = "{:0.0f}".format(float(row[5]))
                row[2] = "match"
                attrib = "ID={tag}:{name_tag}{counter};Name={name}{note}".format(
                    tag=tag,
                    name_tag=name_tag,
                    counter=rep_counter,
                    name=name,
                    note=note,
                )
                print(*row[:8], attrib, sep="\t", file=out_exons, flush=True)
                row[2] = "match_part"
                attrib = "ID={tag}:{name_tag}{counter}-exon1;Parent={tag}:{name_tag}{counter}".format(
                    tag=tag, name_tag=name_tag, counter=rep_counter
                )
                print(*row[:8], attrib, sep="\t", file=out_exons, flush=True)
                print("###", file=out_exons, flush=True)
                row[2] = "exon"
                row[6] = "."
                print(*row[:8], attrib, sep="\t", file=out_exons_unstranded, flush=True)
