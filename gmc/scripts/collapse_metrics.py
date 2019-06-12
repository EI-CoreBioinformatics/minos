import argparse
import sys
import csv

SCORES = [
	"protein_score",
	"transcript_score",
	"hom_qcov_score",
	"hom_tcov_score",
	"hom_acov_score",
	"te_score",
	"cpc_score",
	"expression_score",
	"classification"
]

HEADER = ["#transcript", "gene", "alias"] + SCORES + [s + "_gene" for s in SCORES] + ["confidence", "repeat_associated", "biotype", "discard", "region"]


def initialize_transcript_info(gff):
	transcripts_info = dict()
	
	for row in csv.reader(open(gff), delimiter="\t"):
		if not row[0].startswith("#"):
			if row[2].strip() == "mRNA":
				attrib = dict(item.split("=") for item in row[8].strip(";").split(";"))
				
				if any(map(lambda x: x is None, (attrib.get("ID"), attrib.get("Parent"), attrib.get("Name")))):
					raise ValueError("Error: Cannot parse all variables (ID, Parent, Name). Please check entry:\n{}\n".format("\t".join(row)))

				if not transcripts_info.get(attrib["Name"]) is None:
					raise ValueError("Error: Potential duplicate entry. Transcript '{}' already processed. Please check.\n{}\n".format(attrib["Name"], "\t".join(row)))

				start, end = map(int, row[3:5])
				start, end = (start, end) if start < end else (end, start)

				transcripts_info[attrib["Name"]] = {
					"id": attrib["ID"],
					"parent": attrib["Parent"],
					"alias": attrib["Name"],
					"region": "{}:{}..{}".format(row[0], start, end)
				}

	return transcripts_info

def process_kallisto_data(files):
	def read_kallisto(tsv, kallisto_data):
		for row in csv.DictReader(open(tsv), delimiter="\t"):
			prev_tpm = kallisto_data.get(row["target_id"], 0.0)
			try:
				tpm = float(row["tpm"])
			except:
				tpm = 0.0
			kallisto_data[row["target_id"]] = prev_tpm if prev_tpm > tpm else tpm

	kallisto_data = dict()
	for f in files:
		try:
			read_kallisto(f, kallisto_data)
		except FileNotFoundError:
			print("Warning: Could not find kallisto data file at " + f)

	if not kallisto_data:
		raise ValueError("Error: No kallisto data processed.")

	return kallisto_data

def read_metrics_info(tsv):
	metrics_info = dict()
	for row in csv.reader(open(tsv), delimiter="\t"):
		metric_type = row[0]
		if metric_type == "mikado":
			metric_type = "{}.{}".format(row[0], "protein" if "protein" in row[1] else "transcript")
		metrics_info.setdefault(metric_type, set()).add(row[1])
	return metrics_info


def generate_final_info(metrics_matrix, metrics_info, transcripts_info, kallisto_data):
	model_info, gene_info = dict(), dict()
	for row in csv.DictReader(open(metrics_matrix), delimiter="\t"):

		scores = {
			"protein_score": max(float(row[k]) for k in metrics_info.get("mikado.protein", set()) if k.endswith("_aF1")),
			"transcript_score": max(float(row[k]) for k in metrics_info.get("mikado.transcript", set()) if k.endswith("_aF1")),
			"hom_qcov_score": max(float(row[k]) for k in metrics_info.get("blast", set()) if k.endswith("_qCov"),
			"hom_tcov_score": max(float(row[k]) for k in metrics_info.get("blast", set()) if k.endswith("_tCov"),
			"hom_acov_score": (hom_qcov_score + hom_tcov_score) / 2.0,
			"te_score": 0.0 # !TODO, # we get the highest for the te and when we compute for the gene we take the lowest downstream
			"cpc_score": row["cpc"],
			"expression_score": 0.0,
			"classification": 1 #TODO!
		}


		tinfo = transcripts_info.get(tid)
		if tinfo is not None:

			if not model_info.get(tid) is None:                                                                                                   	
				raise ValueError("Error: Potential duplicate entry. Transcript '{}' already processed. Please check.\n{}\n".format(tid, "\t".join(row)))
			
			kallisto_score = kallisto_data.get(tid)
			if kallisto_score is None:
				raise ValueError("Error: Could not extract tpm data for transcript " + tid)

			scores["expression_score"] = kallisto_score

			model_info["id"] = dict()
			model_info["id"].update(tinfo)
			gid = model_info["id"]["gene"] = tinfo["parent"]
			del model_info["id"]["parent"]
			model_info["id"].update(scores)
			
			# get the highest metrics value for gene
			# except for te_score
			ginfo = gene_info.get(gid, dict())			
			if not ginfo:
				gene_info[gid] = scores
			else:
				for k, v in ginfo:
					cmp_f = max if k != "te_score" else min
					gene_info[gid][k] = cmp_f(v, scores[k])


	return model_info, gene_info


def write_scores(model_info, gene_info):
	
	print(*HEADER, sep="\t")

	for tid, tinfo in sorted(model_info.items(), key=lambda x:x[0]):
		row = [tinfo["id"], tinfo["gene"], tinfo["alias"]]
		row.extend(tinfo.get(score, ".") for score in SCORES)
		gid = tinfo["gene"]
		gene_scores = list(
			gene_info.get(gid, dict()).get(score, ".") for score in SCORES
		)
		if any(score == "." for score in gene_scores):
			raise ValueError("Error: Cannot find all gene scores:\n{}".format("\n".join(zip(SCORES, gene_scores))))
		
		high_confidence = gene_info[gid]["classification"] == 1 or gene_info[gid]["hom_acov_score"] >= 0.8 or (gene_info[gid]["hom_acov_score"] >= 0.6 and gene_info[gid]["transcript_score"] >= 0.4)
		repeat_associated = gene_info[gid]["te_score"] >= 0.4
		biotype = None
		if repeat_associated:
			biotype = "transposable_element_gene"
		elif gene_info[gid]["hom_acov_score"] < 0.3 and gene_info[gid]["cpc_score"] < 0.25 and gene_info[gid]["classification"] == 0:
			biotype = "predicted_gene"
		else:
			biotype = "protein_coding_gene"
		if biotype is None:
			raise ValueError("Error: Could not determine biotype for transcript " + tid)

		discard = not any(gene_info[gid][score] > 0 for score in ("protein_score", "transcript_score", "classification", "hom_acov_score")) and gene_info[gid]["expression_score"] < 0.3
	

		row.extend([
			"High" if high_confidence else "Low",
			str(repeat_associated),
			biotype,
			str(discard),
			tinfo["region"]
		])

		print(*row, sep="\t")


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("input_gff", type=str)
	ap.add_argument("metrics_matrix", type=str)
	ap.add_argument("metrics_info", type=str)
	ap.add_argument("kallisto_tpm", nargs="*")
	args = ap.parse_args()


	print(args)
	try:
		transcripts_info = initialize_transcripts_info(args.input_gff)
	except FileNotFoundError:
		print("Error: Cannot find input gff at " + args.input_gff)

	try:
		metrics_info = read_metrics_info(args.metrics_info)
	except FileNotFoundError:
		print("Error: Cannot find metrics info at " + args.metrics_info)

	kallisto_data = process_kallisto_data(args.kallisto_tpm)

	try:
		model_info, gene_info = generate_final_info(args.metrics_matrix, metrics_info, transcripts_info, kallisto_data)
	except FileNotFoundError:
		print("Error: Cannot find metrics matrix at " + args.metrics_matrix)


	write_scores(model_info, gene_info)


if __name__ == "__main__":
	main()

