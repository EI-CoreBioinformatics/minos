import sys
import os
import argparse
import csv
from collections import Counter

class GffReleaseGenerator():
	def __init__(self, metrics_tsv):
		try:
			self.__read_metrics_info(metrics_tsv)
		except FileNotFoundError:
			raise FileNotFoundError("Could not find metrics data at " + metrics_tsv)

	def __read_metrics_info(self, metrics_tsv):
		transcript_cols = ["#transcript", "gene", "alias", "confidence", "biotype", "discard"]
		gene_cols = ["gene", "confidence", "biotype", "discard"]
	
		self.transcript_metrics_info, self.gene_metrics_info = dict(), dict()
		for row in csv.DictReader(open(metrics_tsv), delimiter="\t"):
			tid, gid = row["#transcript"], row["gene"]
			if self.transcript_metrics_info.get(tid, None) is not None:
				raise ValueError("Error: Transcript '{}' has already been processed - caused by potential duplicate entry.".format(tid))
	
			self.transcript_metrics_info[tid] = {
				col: row[col] for col in transcript_cols
			}
		
			if self.gene_metrics_info.get(gid, None) is None:
				self.gene_metrics_info[gid] = {
					col: row[col] for col in gene_cols
				}

	
	gff_header = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
	valid_feature_types = {"gene", "mRNA", "exon", "CDS", "five_prime_UTR", "three_prime_UTR"}
	
	def process_gff(self, raw_gff, name_prefix="XYZ_v1", file_prefix="mikado.annotation"):

		gene_counter = 0
		new_genes = dict()
		new_transcripts = dict()

		def get_attrib(key, attributes, feature):
			try:
				return attributes[key]
			except:
				raise ValueError("Fatal Error: Cannot parse {} {} field.".format(feature, key))

		out_release = os.path.join(os.path.dirname(raw_gff), file_prefix + ".release.unsorted.gff3")
		out_browser = os.path.join(os.path.dirname(raw_gff), file_prefix + ".release_browser.unsorted.gff3")

		with open(raw_gff + ".old_new_id_relation.txt", "wt") as logfile, open(out_release, "wt") as out_rel, open(out_browser, "wt") as out_brw:
			print(*["#New_Gene_ID", "New_mRNA_ID", "Old_Gene_ID", "Old_mRNA_ID"], sep="\t", file=logfile)
			feature_counter = Counter()

			try:			
				for row in csv.DictReader(open(raw_gff), delimiter="\t", fieldnames=self.gff_header):
					if row["seqid"].startswith("#"):
						print(*(c for c in row.values() if c is not None), sep="\t", file=out_rel)
						print(*(c for c in row.values() if c is not None), sep="\t", file=out_brw)
					else:
						start, end = int(row["start"]), int(row["end"])
						start, end = (start, end) if start < end else (end, start)
						attrib = dict(item.split("=") for item in row["attributes"].strip(";").split(";"))

						row["source"] = name_prefix

						if row["type"] == "gene":
							gid = get_attrib("ID", attrib, row["type"])

							ginfo = self.gene_metrics_info.get(gid, None)
							if ginfo is None:
								raise ValueError("Error: Gene '{}' is not in the metrics file.\n{}\n".format(gid, "\t".join(row.values())))
					
							if eval(ginfo["discard"]):
								print("Debug: Gene discarded '{}'".format(gid), file=sys.stderr)
								continue

							gene_counter += 10
							new_gene_id = "{}_{:07d}".format(name_prefix, gene_counter)
							
							if new_genes.get(gid, None) is not None:
								raise ValueError("Error: Duplicated gene id enocuntered in the input gff '{}'".format(gid))
							new_genes[gid] = new_gene_id
							
							row["attributes"] = "ID={0};Name={0};biotype={1};confidence={2}".format(new_gene_id, ginfo["biotype"], ginfo["confidence"])

							print(*row.values(), sep="\t", file=out_rel)
							print(*row.values(), sep="\t", file=out_brw)

						elif row["type"] == "mRNA":
							tid = get_attrib("ID", attrib, row["type"])
							gid = get_attrib("Parent", attrib, row["type"])

							# note = eval(get_attrib("Note", attrib, row["type"]))
							note = get_attrib("Note", attrib, row["type"]).split(",")
							note = dict(item.split(":") for item in note if item.count(":") == 1)
							is_primary = note.get("primary") is not None and eval(note.get("primary"))
							region = "{}:{}..{}".format(row["seqid"], start, end)

							tinfo = self.transcript_metrics_info.get(tid, None)
							if tinfo is None:
								raise ValueError("Error: Transcript '{}' is not in the metrics file.\n{}\n".format(tid, "\t".join(row.values())))

							if eval(tinfo["discard"]):
								print("Debug: Transcript discarded '{}'".format(tid), file=sys.stderr)
								continue

							try:
								# Gemy's original:
								# my ($count_mRNA) = $mrna =~ /\S+\D\S+\D+(\d+)/; 
								# get the last 1 from mikado.109676G2.1, I might need to change this for different annotation
								primary_suffix_check = tid.split(".")[-1]
							except:
								raise ValueError("Error: Cannot extract suffix count from transcript id ({}).".format(tid))

							new_gene = new_genes.get(gid, None)
							if new_gene is None:
								raise ValueError("Error: No gene feature for transcript id {} (gene={}).".format(tid, gid))

							new_transcript = "{}.{}".format(new_gene, primary_suffix_check)
							if new_transcripts.get(tid, None) is not None:
								raise ValueError("Error: Duplicated transcript id '{}' ({})".format(tid, new_transcript))

							new_transcripts[tid] = new_transcript
							feature_counter = Counter()

							attribs = "ID={0};Parent={1};Name={0};Note={1}".format(new_transcript, new_gene)
							# browser gff has different attributes
							row["attributes"] = attribs + "|{}|conf:{}|rep:{}".format(tinfo["biotype"], tinfo["confidence"], is_primary)
							print(*row.values(), sep="\t", file=out_brw)
							# release gff
							row["attributes"] = attribs + ";confidence={};representative={}".format(tinfo["confidence"], is_primary)
							print(*row.values(), sep="\t", file=out_rel)

							print(new_gene, new_transcript, gid, tid, sep="\t", file=logfile)

						elif row["type"] in self.valid_feature_types:
						
							feature_counter[row["type"]] += 1

							tid = get_attrib("Parent", attrib, row["type"])
							tinfo = self.transcript_metrics_info.get(tid, None)
							if tinfo is None:
								raise ValueError("Error: Transcript '{}' is not in the metrics file.\n{}\n".format(tid, "\t".join(row.values())))

							if eval(tinfo["discard"]):
								print("Debug: {} parent transcript discarded '{}'".format(row["type"], tid), file=sys.stderr)
								continue

							new_transcript = new_transcripts.get(tid, None)
							if new_transcript is not None:
								new_feature_id = "{}.{}{}".format(new_transcript, row["type"], feature_counter[row["type"]])
								row["attributes"] = "ID={};Parent={}".format(new_feature_id, new_transcript)
								print(*row.values(), sep="\t", file=out_rel)
								print(*row.values(), sep="\t", file=out_brw)

						else:
							print("WARN: Unknown feature type '{}'\n{}\n".format(row["type"], "\t".join(row.values())), file=sys.stderr)

			except FileNotFoundError:
				raise FileNotFoundError("Error: Could not find gff file at {}".format(raw_gff))


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("input_gff", type=str)
	ap.add_argument("collapsed_metrics", type=str)
	ap.add_argument("--annotation-version", type=str, default="EIv1")
	ap.add_argument("--genus-identifier", type=str, default="XYZ")
	args = ap.parse_args()


	gff_rg = GffReleaseGenerator(args.collapsed_metrics).process_gff(
		args.input_gff, 
		name_prefix="{}_{}".format(args.genus_identifier, args.annotation_version)
	)
	
	

	
	

if __name__ == "__main__": 
	main()
