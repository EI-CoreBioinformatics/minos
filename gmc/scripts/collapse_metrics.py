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

# 3 + 9 * 2 + 5 = 26
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
		metric_type, metric_id = row[0:2]
		metrics_info.setdefault(metric_type, set()).add(metric_id)
	return metrics_info


def generate_final_info(metrics_matrix, metrics_info, transcripts_info, kallisto_data):
	model_info, gene_info = dict(), dict()
	for row in csv.DictReader(open(metrics_matrix), delimiter="\t"):
		tid = row["tid"]

		scores = {
			"protein_score": max(float(row[k + "_aF1"]) for k in metrics_info.get("mikado.protein", set())),
			"transcript_score": max(float(row[k + "_aF1"]) for k in metrics_info.get("mikado.transcript", set())),
			"hom_qcov_score": max(float(row[k + "_qCov"]) for k in metrics_info.get("blast", set())),
			"hom_tcov_score": max(float(row[k + "_tCov"]) for k in metrics_info.get("blast", set())),
			"hom_acov_score": 0, 
			"te_score": 0.0,  # !TODO, # we get the highest for the te and when we compute for the gene we take the lowest downstream
			"cpc_score": float(row["cpc"]),
			"expression_score": 0.0,
			"classification": 0 ## cschu 20200203: issue9: disabled until full-lengther replacement implemented
		}
		scores["hom_acov_score"] = (scores["hom_qcov_score"] + scores["hom_tcov_score"]) / 2.0

		tinfo = transcripts_info.get(tid, None)
		if tinfo is not None:
			# print(tid, tinfo, file=sys.stderr)

			if not model_info.get(tid) is None:                                                                                                   	
				raise ValueError("Error: Potential duplicate entry. Transcript '{}' already processed. Please check.\n{}\n".format(tid, "\t".join(row)))
			
			kallisto_score = kallisto_data.get(tinfo["id"])
			if kallisto_score is None:
				raise ValueError("Error: Could not extract tpm data for transcript {} ({})".format(tid, tinfo["id"]))

			scores["expression_score"] = kallisto_score

			model_info[tid] = dict()
			model_info[tid].update(tinfo)
			gid = model_info[tid]["gene"] = tinfo["parent"]
			del model_info[tid]["parent"]
			model_info[tid].update(scores)
			
			# get the highest metrics value for gene
			# except for te_score
			ginfo = gene_info.get(gid, dict())			
			if not ginfo:
				gene_info[gid] = scores
			else:
				for k, v in ginfo.items():
					cmp_f = max if k != "te_score" else min
					gene_info[gid][k] = cmp_f(v, scores[k])


	return model_info, gene_info


def cmp_score(a, b, op):
	if op == "eq":
		return a == b
	if op == "lt":
		return a < b
	if op == "gt":
		return a > b
	if op == "ge":
		return a >= b
	raise ValueError("Invalid check: {} {} {}".format(a, b, op))

def check_expression(expression, values):
	if type(expression[0]) is str:
		a, b, op = expression
		expression = (values[a], b, op)
		return cmp_score(*expression)
	operator = expression[0]
	return operator(
		check_expression(exp, values) for exp in expression[1:]
	)

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
		row.extend(gene_scores)
		
		biotype = "protein_coding_gene"
		repeat_associated = gene_info[gid]["te_score"] >= 0.4
		if repeat_associated:
			biotype = "transposable_element_gene"
		else:
			predicted_gene_checks = (
				all,
				("hom_acov_score", 0.3, "lt"), 
				("cpc_score", 0.25, "lt"),
				# and gene_info[gid]["classification"] == 0: ## cschu 20200203: issue9: disabled until full-lengther replacement implemented
			)
			if check_expression(predicted_gene_checks, gene_info[gid]):
				biotype = "predicted_gene"

		hiconf_checks = (
			any,
			("classification", 1, "eq"),
			("hom_acov_score", 0.8, "ge"),
			(all, ("hom_acov_score", 0.6, "ge"), ("transcript_score", 0.4, "ge"))
		)
		high_confidence = check_expression(hiconf_checks, gene_info[gid])

		discard_checks = (
			all,
			("protein_score", 0, "eq"),
			("transcript_score", 0, "eq"),
			("hom_acov_score", 0, "eq"),
			("expression_score", 0.3, "lt"),
			# ("classification", "eq", 0) ## cschu 20200203: issue9: disabled until full-lengther replacement implemented 
		)
		discard = check_expression(discard_checks, gene_info[gid])

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


	# print(args)
	try:
		transcripts_info = initialize_transcript_info(args.input_gff)
	except FileNotFoundError:
		print("Error: Cannot find input gff at " + args.input_gff, file=sys.stderr)

	try:
		metrics_info = read_metrics_info(args.metrics_info)
	except FileNotFoundError:
		print("Error: Cannot find metrics info at " + args.metrics_info, file=sys.stderr)

	kallisto_data = process_kallisto_data(args.kallisto_tpm)

	try:
		model_info, gene_info = generate_final_info(args.metrics_matrix, metrics_info, transcripts_info, kallisto_data)
	except FileNotFoundError:
		print("Error: Cannot find metrics matrix at " + args.metrics_matrix, file=sys.stderr)


	write_scores(model_info, gene_info)


if __name__ == "__main__":
	main()

