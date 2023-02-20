import csv
import re


def extract_exons(_in, _out):
    with open(_out, "w") as exons_out:
        exon = 1
        regex = re.compile('([^ ;]+) +("?[^\'"]+"?) *;?')
        for row in csv.reader(open(_in), delimiter="\t"):
            if row and not row[0].startswith("#") and row[2].lower() == "exon":
                attr = dict(
                    (item.group(1).strip(), item.group(2).strip())
                    for item in regex.finditer(row[8])
                )
                row[8] = 'ID="{tid}.exon{exon}";Parent="{tid}";'.format(
                    tid=attr["transcript_id"].strip('"'), exon=exon
                )
                exon += 1
                print(*row, sep="\t", flush=True, file=exons_out)


#        """
#        awk '$3 == "exon"' {input[0]} | awk -v FS="\\t" -v OFS="\\t" '{{exon += 1; split($9,a,"\\""); $9="ID="a[4]".exon"exon";Parent="a[4]; print $0}}' > {output[0]}
#        """.strip().replace("\n\t", " ")
