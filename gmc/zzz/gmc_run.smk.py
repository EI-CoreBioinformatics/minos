import os
import sys
from enum import Enum, unique, auto

EXTERNAL_METRICS_DIR = os.path.join(config["outdir"], "generate_metrics")
LOG_DIR = os.path.join(config["outdir"], "logs")
TEMP_DIR = os.path.join(config["outdir"], "tmp")

@unique
class ExternalMetrics(Enum):
	MIKADO_TRANSCRIPTS_OR_PROTEINS = auto()
	PROTEIN_BLAST_TOPHITS = auto()
	CPC_CODING_POTENTIAL = auto()
	KALLISTO_TPM_EXPRESSION	= auto()
	REPEAT_ANNOTATION = auto()



def get_rnaseq(wc):
	return [item for sublist in config["data"]["expression-runs"][wc.run] for item in sublist]

def get_all_transcript_assemblies(wc):
	print([config["data"]["transcript-runs"][asm][0] for asm in config["data"]["transcript-runs"]])
	return [config["data"]["transcript-runs"][asm][0] for asm in config["data"]["transcript-runs"]]

def get_protein_alignments(wc):
	return config["data"]["protein-runs"][wc.run][0]
def get_transcript_alignments(wc):
	return config["data"]["transcript-runs"][wc.run][0]
def get_protein_sequences(wc):
	return config["data"]["protein-seqs"].get(wc.run, [""])[0]

OUTPUTS = [
	os.path.join(config["outdir"], "mikado_prepared.fasta"), 
	os.path.join(EXTERNAL_METRICS_DIR, "CPC-2.0_beta", "mikado_prepared.fasta.cpc2output.txt"),
]

if config["use-tpm-for-picking"]:
	OUTPUTS.append(
		os.path.join(EXTERNAL_METRICS_DIR, "kallisto", "mikado_prepared.fasta.idx")
	)

POST_PICK_EXPRESSION = list()
for run in config.get("data", dict()).get("expression-runs", dict()):
	if config["use-tpm-for-picking"]:
		OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "kallisto", run, "abundance.tsv"))
	POST_PICK_EXPRESSION.append(os.path.join(config["outdir"], "kallisto", run, "abundance.tsv"))

# gmc_run4/generate_metrics/mikado_compare/vs_transcripts/Scallop_Old_leaf_transcripts/vs_transcripts_Scallop_Old_leaf_transcripts.refmap
for run in config.get("data", dict()).get("protein-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "proteins", run, "mikado_" + run + ".refmap"))
for run in config.get("data", dict()).get("transcript-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "transcripts", run, "mikado_" + run + ".refmap"))
for run in config.get("data", dict()).get("protein-seqs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], run, run + ".{}.tsv.tophit".format(config["blast-mode"])))


localrules:
	all,
	gmc_metrics_blastp_combine,
	gmc_metrics_generate_metrics_info,
	gmc_metrics_generate_metrics_matrix,
	gmc_parse_mikado_pick,
	gmc_gffread_extract_proteins_post_pick,
	gmc_gffread_extract_proteins,
	gmc_protein_completeness,
	gmc_gffread_extract_cdna_post_pick,
	gmc_gff_genometools_check_post_pick,
	gmc_collect_biotype_conf_stats

rule all:
	input:
		OUTPUTS,
		os.path.join(EXTERNAL_METRICS_DIR, "metrics_info.txt"),
		os.path.join(config["outdir"], "MIKADO_SERIALISE_DONE"),
		os.path.join(config["outdir"], "mikado.subloci.gff3"),
		os.path.join(config["outdir"], "mikado.loci.gff3"),
		os.path.join(config["outdir"], "mikado.annotation.gff"),
		os.path.join(config["outdir"], "mikado.annotation.proteins.fasta"),
		os.path.join(config["outdir"], "mikado.annotation.protein_status.tsv"),
		os.path.join(config["outdir"], "mikado.annotation.protein_status.summary"),
		os.path.join(config["outdir"], "mikado.annotation.cdna.fasta"),
		os.path.join(config["outdir"], "kallisto", "mikado.annotation.cdna.fasta.idx"),
		POST_PICK_EXPRESSION,
		os.path.join(config["outdir"], "mikado.annotation.gt_checked.gff"),
		os.path.join(config["outdir"], "mikado.annotation.collapsed_metrics.tsv"),
		os.path.join(config["outdir"], "mikado.annotation.gt_checked.validation_report.txt"),
		os.path.join(config["outdir"], "mikado.annotation.release.unsorted.gff3"),
		os.path.join(config["outdir"], "mikado.annotation.release_browser.unsorted.gff3"),
		os.path.join(config["outdir"], config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1") + ".release.gff3"),
		os.path.join(config["outdir"], config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1") + ".release_browser.gff3"),  
		os.path.join(config["outdir"], config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1") + ".sanity_checked.release.gff3"),
		os.path.join(config["outdir"], config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1") + ".sanity_checked.release.gff3.mikado_stats.txt"),
		os.path.join(config["outdir"], config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1") + ".sanity_checked.release.gff3.mikado_stats.tsv"),
		os.path.join(config["outdir"], config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1") + ".sanity_checked.release.gff3.biotype_conf.tsv"),
		[
			os.path.join(config["outdir"], config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1") + ".sanity_checked.release.gff3") + ".{}.fasta".format(dtype) for dtype in {"cdna", "cds", "pep.raw", "pep"}
		],
		[os.path.join(config["outdir"], config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1") + ".sanity_checked.release.gff3") + ".pep.raw.{}".format(suffix) for suffix in {"protein_status.tsv", "protein_status.summary"}]	

rule gmc_mikado_prepare:
	input:
		config["mikado-config-file"]
	output:
		os.path.join(config["outdir"], "mikado_prepared.fasta"),
		os.path.join(config["outdir"], "mikado_prepared.gtf")
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="prepare"),
		program_params = config["params"]["mikado"]["prepare"],
		outdir = config["outdir"]
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_prepare.log")
	threads:
		30
	shell:
		"{params.program_call} {params.program_params} --json-conf {input[0]} --procs {threads} -od {params.outdir} &> {log}"

rule gmc_mikado_compare_index_reference:
	input:
		rules.gmc_mikado_prepare.output[1]
	output:
		rules.gmc_mikado_prepare.output[1] + ".midx"
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="compare"),
		program_params = config["params"]["mikado"]["compare"]["index"]
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_compare_index_reference.log")
	shell:
		"{params.program_call} {params.program_params} -r {input} &> {log}"

rule gmc_gffread_extract_proteins:
	input:
		gtf = rules.gmc_mikado_prepare.output[1],
		refseq = config["reference-sequence"]
	output:
		rules.gmc_mikado_prepare.output[1] + (".prot.fasta" if config["blast-mode"] == "blastp" else ".cds.fasta")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".gffread_extract.log")
	threads:
		1
	params:
		program_call = config["program_calls"]["gffread"],
		program_params = config["params"]["gffread"][config["blast-mode"]]
	shell:
		"{params.program_call} {input.gtf} -g {input.refseq} {params.program_params} {output[0]}.raw &> {log} && " + \
		"awk '/^[^>]/ {{ $1=gensub(\"\\\\.\", \"\", \"g\", $1) }} {{ print $0 }}' {output[0]}.raw > {output[0]}"


rule gmc_metrics_cpc2:
	input:
		rules.gmc_mikado_prepare.output[0]
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "CPC-2.0_beta", os.path.basename(rules.gmc_mikado_prepare.output[0]) + ".cpc2output.txt")
	params:
		program_call = config["program_calls"]["cpc2"],
		program_params = config["params"]["cpc2"]
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".CPC2.log")
	threads:
		4
	shell:
		"{params.program_call} {params.program_params} -i {input[0]} -o {output[0]} &> {log}"


if config["use-tpm-for-picking"]:

	rule gmc_metrics_kallisto_index:
		input:
			rules.gmc_mikado_prepare.output[0]
		output:
			os.path.join(EXTERNAL_METRICS_DIR, "kallisto", os.path.basename(rules.gmc_mikado_prepare.output[0]) + ".idx")
		log:
			os.path.join(LOG_DIR, config["prefix"] + ".kallisto_index.log")
		params:
			program_call = config["program_calls"]["kallisto"].format(program="index"),
			program_params = config["params"].get("kallisto", {}).get("index", "")
		shell:
			"{params.program_call} {params.program_params} -i {output[0]} {input[0]} &> {log}"

	rule gmc_metrics_kallisto_quant:
		input:
			index = rules.gmc_metrics_kallisto_index.output[0],
			reads = get_rnaseq
		output:
			os.path.join(EXTERNAL_METRICS_DIR, "kallisto", "{run}", "abundance.tsv")
		log:
			os.path.join(LOG_DIR, config["prefix"] + ".{run}.kallisto.log")
		params:
			program_call = config["program_calls"]["kallisto"].format(program="quant"),
			program_params = config["params"].get("kallisto", {}).get("quant", ""),
			stranded = lambda wildcards: "" if wildcards.run.endswith("_xx") else "--" + wildcards.run[-2:] + "-stranded",
			outdir = os.path.join(EXTERNAL_METRICS_DIR, "kallisto", "{run}")
		threads:
			32
		shell:
			"{params.program_call} {params.program_params} {params.stranded} -i {input.index} -o {params.outdir} --threads {threads} {input.reads} &> {log}"

rule gmc_metrics_mikado_compare_vs_transcripts:
	input:
		midx = rules.gmc_mikado_compare_index_reference.output[0],
		mika = rules.gmc_mikado_prepare.output[1],
		transcripts = get_transcript_alignments
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "transcripts", "{run}", "mikado_{run}.refmap")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_compare.tran.{run}.log")
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="compare"),
		program_params = config["params"]["mikado"]["compare"]["transcripts"],
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "transcripts", wildcards.run),
		transcripts = lambda wildcards: wildcards.run
	shell:
		"mkdir -p {params.outdir} && " + \
		"{params.program_call} {params.program_params} -r {input.mika} -p {input.transcripts} -o {params.outdir}/mikado_{params.transcripts} &> {log} && " + \
		"touch {output[0]}"
		
rule gmc_metrics_mikado_compare_vs_proteins:
	input:
		midx = rules.gmc_mikado_compare_index_reference.output[0],
		mika = rules.gmc_mikado_prepare.output[1],
		proteins = get_protein_alignments
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "proteins", "{run}", "mikado_{run}.refmap")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_compare.prot.{run}.log")
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="compare"),
		program_params = config["params"]["mikado"]["compare"]["proteins"],
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "proteins", wildcards.run),
		proteins = lambda wildcards: wildcards.run
	shell:
		"mkdir -p {params.outdir} && " + \
		"{params.program_call} {params.program_params} -r {input.mika} -p {input.proteins} -o {params.outdir}/mikado_{params.proteins} &> {log} && " + \
		"touch {output[0]}"

rule gmc_metrics_blastp_mkdb:
	input:
		get_protein_sequences
	output:
		os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], "{run}", "blastdb", "{run}")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".makeblastdb.{run}.log")
	params:
		program_call = config["program_calls"]["blast"].format(program="makeblastdb"),
		program_params = config["params"]["blast"]["makeblastdb"], 
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], wildcards.run, "blastdb"),
		db_prefix = lambda wildcards: wildcards.run
	shell:
		"{params.program_call} {params.program_params} -in {input[0]} -out {params.outdir}/{params.db_prefix} -logfile {log} && " + \
		"touch {output[0]}"

checkpoint gmc_chunk_proteins:
	input:
		rules.gmc_gffread_extract_proteins.output[0]
	output:
		chunk_dir = directory(TEMP_DIR)
	log:
		os.path.join(LOG_DIR, os.path.basename(rules.gmc_gffread_extract_proteins.output[0]) + ".chunk.log")
	params:
		chunksize = 1000,
		outdir = TEMP_DIR
	shell:
		"mkdir -p {params.outdir} && " + \
		# awk script by Pierre Lindenbaum https://www.biostars.org/p/13270/
		"awk 'BEGIN {{n=0;m=1;}} /^>/ {{ if (n%{params.chunksize}==0) {{f=sprintf(\"{params.outdir}/chunk-%d.txt\",m); m++;}}; n++; }} {{ print >> f }}' {input[0]} &> {log}"		

rule gmc_metrics_blastp_chunked:
	input:
		chunk = os.path.join(TEMP_DIR, "chunk-{chunk}.txt"),
		db = rules.gmc_metrics_blastp_mkdb.output[0]
	output:
		os.path.join(TEMP_DIR, "{run}", "chunk-{chunk}." + config["blast-mode"] + ".tsv")
	log:
		os.path.join(LOG_DIR, "blast_logs", "chunk-{chunk}.{run}." + config["blast-mode"] + ".log")
	threads:
		8
	params:
		program_call = config["program_calls"]["blast"].format(program=config["blast-mode"]),
		program_params = config["params"]["blast"][config["blast-mode"]]
	shell:
		"{params.program_call} {params.program_params} -query {input.chunk} -out {output[0]} -num_threads {threads} " + \
		"-db {input.db} -outfmt \"6 qseqid sseqid pident qstart qend sstart send qlen slen length nident mismatch positive gapopen gaps evalue bitscore\" &> {log}"


def aggregate_input(wildcards):
	checkpoint_output = checkpoints.gmc_chunk_proteins.get(**wildcards).output.chunk_dir
	return expand(
		os.path.join(TEMP_DIR, "{run}", "chunk-{chunk}." + config["blast-mode"] + ".tsv"),
		run=wildcards.run,
		chunk=glob_wildcards(os.path.join(checkpoint_output, "chunk-{chunk}.txt")).chunk
	)


rule gmc_metrics_blastp_combine:
	input:
		aggregate_input 
	output:
		os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], "{run}", "{run}." + config["blast-mode"] + ".tsv")		
	shell:
		"cat {input} > {output[0]} && rm {input}"


rule gmc_metrics_blastp_tophit:
	input:
		rules.gmc_metrics_blastp_combine.output[0]
	output:
		rules.gmc_metrics_blastp_combine.output[0] + ".tophit"
	params:
		pident_threshold = config["params"]["blast"]["tophit"]["pident_threshold"],
		qcov_threshold = config["params"]["blast"]["tophit"]["qcov_threshold"]
	run:
		import csv
		def calc_coverage_perc(start, end, length):
			alen = (end - start + 1) if start < end else (start - end + 1)
			return round(alen/length * 100, 2)
		with open(output[0], "w") as blast_out:
			seen = set()
			for row in csv.reader(open(input[0]), delimiter="\t"):
				if not row or row[0].startswith("#") or row[0] in seen:
					continue
				if len(row) < 17:
					raise ValueError("Misformatted blastp line (ncols={ncols}) in {blastp_input}. Expected is outfmt 6 (17 cols).".format(
						ncols=len(row), 
						blastp_input=input[0])
					)
				try:
					qstart, qend, sstart, send, qlen, slen = map(int, row[3:9])
				except:
					raise ValueError("Misformatted blastp line: could not parse integer values from columns 2-7 ({})".format(",".join(row[3:9])))
				qcov, scov = calc_coverage_perc(qstart, qend, qlen), calc_coverage_perc(sstart, send, slen)
				if float(row[2]) >= params.pident_threshold and qcov >= params.qcov_threshold:
					print(*row, "{:.2f}".format(qcov), "{:.2f}".format(scov), sep="\t", file=blast_out)
					seen.add(row[0])


rule gmc_metrics_generate_metrics_info:
	input:
		prev_outputs = OUTPUTS
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "metrics_info.txt")
	run:
		import os
		import glob
		import csv

		# this block is unstable against tempering with the structure of the generate_metrics dir
		# might work better as state-machine
		with open(output[0], "wt") as metrics_info:
			walk = os.walk(EXTERNAL_METRICS_DIR)
			next(walk)
			rows = list()

			while walk:
				try:
					cwd, dirs, files = next(walk)
				except StopIteration:
					break
				cwd_base = os.path.basename(cwd)
				if cwd_base == "CPC-2.0_beta":
					mclass, mid, path = "cpc", "cpc", os.path.join(cwd, files[0])
					rows.append((ExternalMetrics.CPC_CODING_POTENTIAL, mclass, mid, path))
				elif cwd_base in {"proteins", "transcripts"}:
					mclass = "mikado.{}".format(cwd_base[:-1])
					for mid in dirs:
						cwd, _, files = next(walk)
						path = glob.glob(os.path.join(cwd, "*.refmap"))[0]
						rows.append((ExternalMetrics.MIKADO_TRANSCRIPTS_OR_PROTEINS, mclass, mid, path))
				elif cwd_base in {"blastp", "blastx"}:
					mclass = "blast"
					for mid in dirs:
						cwd, _, files = next(walk)
						path = glob.glob(os.path.join(cwd, "*.tophit"))[0]
						rows.append((ExternalMetrics.PROTEIN_BLAST_TOPHITS, mclass, mid, path))
						# last block is not necessary if we clean up the blast databases before
						try:
							_ = next(walk)
						except StopIteration:
							pass
				elif cwd_base == "kallisto":
					mclass = "expression"
					for mid in dirs:
						cwd, _, _ = next(walk)
						path = os.path.join(cwd, "abundance.tsv")
						rows.append((ExternalMetrics.KALLISTO_TPM_EXPRESSION, mclass, mid, path))

			for mid in config["data"].get("repeat-data", dict()):
				mclass = "repeat"
				path = config["data"]["repeat-data"][mid][0][0]
				rows.append((ExternalMetrics.REPEAT_ANNOTATION, mclass, mid, path))

			for _, mclass, mid, path in sorted(rows, key=lambda x:x[0].value):
			
				print(mclass, mid, path, sep="\t", file=sys.stderr)
				print(mclass, mid, os.path.abspath(path), sep="\t", file=metrics_info)

rule gmc_metrics_generate_metrics_matrix:
	input:
		rules.gmc_metrics_generate_metrics_info.output[0]
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "metrics_matrix.txt")
	log:
		os.path.join(LOG_DIR, "generate_metrics_matrix.log")
	shell:
		"generate_metrics {input[0]} > {output[0]} 2> {log}"

rule gmc_mikado_serialise:
	input:
		config = config["mikado-config-file"],
		ext_scores = rules.gmc_metrics_generate_metrics_matrix.output[0],
		transcripts = rules.gmc_mikado_prepare.output[0]
	output:
		os.path.join(config["outdir"], "MIKADO_SERIALISE_DONE"),
		os.path.join(config["outdir"], "mikado.db")
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="serialise"),
		program_params = config["params"]["mikado"]["serialise"],
		outdir = config["outdir"]
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_serialise.log")
	threads:
		32
	shell:
		"{params.program_call} {params.program_params} --transcripts {input.transcripts} --external-scores {input.ext_scores} --json-conf {input.config} --procs {threads} -od {params.outdir} &> {log} && " + \
		"touch {output[0]}"

rule gmc_mikado_pick:
	input:
		config = config["mikado-config-file"],
		gtf = rules.gmc_mikado_prepare.output[1],
		serialise_done = rules.gmc_mikado_serialise.output[0],
		db = rules.gmc_mikado_serialise.output[1]
	output:
		loci = os.path.join(config["outdir"], "mikado.loci.gff3"),
		subloci = os.path.join(config["outdir"], "mikado.subloci.gff3")
	threads:
		30
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="pick"),
		program_params = config["params"]["mikado"]["pick"],
		outdir = config["outdir"]
	shell:
		"{params.program_call} {params.program_params} -od {params.outdir} --procs {threads} --json-conf {input.config} --subloci-out $(basename {output.subloci}) -db {input.db} {input.gtf}"

rule gmc_parse_mikado_pick:
	input:
		loci = rules.gmc_mikado_pick.output[0]
	output:
		gff = os.path.join(config["outdir"], "mikado.annotation.gff")
	shell:
		"parse_mikado_gff {input.loci} > {output.gff}"

rule gmc_gffread_extract_proteins_post_pick:
	input:
		gff = rules.gmc_parse_mikado_pick.output[0],
		refseq = config["reference-sequence"]
	output:
		os.path.join(config["outdir"], "mikado.annotation.proteins.fasta")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".gffread_extract_post_pick.log")
	params:
		program_call = config["program_calls"]["gffread"],
		program_params = config["params"]["gffread"]["default"]
	shell:
		"{params.program_call} {input.gff} -g {input.refseq} {params.program_params} -x {output[0]} &> {log}" 

rule gmc_gffread_extract_cdna_post_pick:
	input:
		gff = rules.gmc_parse_mikado_pick.output[0],
		refseq = config["reference-sequence"]
	output:
		os.path.join(config["outdir"], "mikado.annotation.cdna.fasta")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".gffread_extract_cdna_post_pick.log")
	params:
		program_call = config["program_calls"]["gffread"],
		program_params = config["params"]["gffread"]["default"]
	shell:
		"{params.program_call} {input.gff} -g {input.refseq} {params.program_params} -w {output[0]} &> {log}" 

rule gmc_gff_genometools_check_post_pick:
	input:
		rules.gmc_parse_mikado_pick.output[0]
	output:
		gff = os.path.join(config["outdir"], "mikado.annotation.gt_checked.gff")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".post_pick_genometools_check.log")
	params:
		program_call = config["program_calls"]["genometools"],
		program_params = config["params"]["genometools"]["check"]
	shell:
		"{params.program_call} {params.program_params} {input[0]} > {output.gff} 2> {log}"

rule gmc_gff_validate_post_gt:
	input:
		rules.gmc_gff_genometools_check_post_pick.output[0]
	output:
		os.path.join(config["outdir"], "mikado.annotation.gt_checked.validation_report.txt")
	shell:
		"validate_gff3 {input} > {output}"

rule gmc_kallisto_index_post_pick:
	input:
		rules.gmc_gffread_extract_cdna_post_pick.output[0]
	output:
		os.path.join(config["outdir"], "kallisto", os.path.basename(rules.gmc_gffread_extract_cdna_post_pick.output[0]) + ".idx")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".kallisto_index_post_pick.log")
	params:
		program_call = config["program_calls"]["kallisto"].format(program="index"),
		program_params = config["params"].get("kallisto", {}).get("index", "")
	shell:
		"{params.program_call} {params.program_params} -i {output[0]} {input[0]} &> {log}"

rule gmc_kallisto_quant_post_pick:
	input:
		index = rules.gmc_kallisto_index_post_pick.output[0],
		reads = get_rnaseq
	output:
		os.path.join(config["outdir"], "kallisto", "{run}", "abundance.tsv")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".{run}.kallisto.post_pick.log")
	params:
		program_call = config["program_calls"]["kallisto"].format(program="quant"),
		program_params = config["params"].get("kallisto", {}).get("quant", ""),
		stranded = lambda wildcards: "" if wildcards.run.endswith("_xx") else "--" + wildcards.run[-2:] + "-stranded",
		outdir = os.path.join(config["outdir"], "kallisto", "{run}")
	threads:
		32
	shell:
		"{params.program_call} {params.program_params} {params.stranded} -i {input.index} -o {params.outdir} --threads {threads} {input.reads} &> {log}"

rule gmc_collapse_metrics:
	input:
		gff = rules.gmc_parse_mikado_pick.output[0],
		ext_scores = rules.gmc_metrics_generate_metrics_matrix.output[0],
		metrics_info = rules.gmc_metrics_generate_metrics_info.output[0],
		expression = expand(rules.gmc_kallisto_quant_post_pick.output, run=config["data"]["expression-runs"].keys())
	output:
		os.path.join(config["outdir"], "mikado.annotation.collapsed_metrics.tsv")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".collapse_metrics.log")
	run:
		from gmc.scripts.collapse_metrics import MetricCollapser
		mc = MetricCollapser(input.gff, input.metrics_info, input.ext_scores, input.expression)
		with open(output[0], "w") as out:
			mc.write_scores(config["collapse_metrics_thresholds"], stream=out)

rule gmc_create_release_gffs:
	input:
		gff = rules.gmc_gff_genometools_check_post_pick.output[0],
		metrics_info = rules.gmc_collapse_metrics.output[0]
	output:
		os.path.join(config["outdir"], "mikado.annotation.release.unsorted.gff3"),
		os.path.join(config["outdir"], "mikado.annotation.release_browser.unsorted.gff3")
	params:
		annotation_version = config.get("annotation_version", "EIv1"),
		genus_identifier = config.get("genus_identifier", "XYZ")
	shell:
		"create_release_gff3 {input.gff} {input.metrics_info} --annotation-version {params.annotation_version} --genus-identifier {params.genus_identifier} 2> {LOG_DIR}/create_release_gff.log"

rule gmc_sort_release_gffs:
	input:
		rules.gmc_create_release_gffs.output
	output:
		os.path.join(config["outdir"], config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1") + ".release.gff3"),
		os.path.join(config["outdir"], config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1") + ".release_browser.gff3")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".sort_release_gffs.log")
	params:
		program_call = config["program_calls"]["genometools"],
		program_params = config["params"]["genometools"]["sort"]
	shell:
		"{params.program_call} {params.program_params} {input[0]} > {output[0]} 2> {log} && " + \
		"{params.program_call} {params.program_params} {input[1]} > {output[1]} 2>> {log}" 

rule gmc_final_sanity_check:
	input:
		rules.gmc_sort_release_gffs.output[0]
	output:
		os.path.join(config["outdir"], config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1") + ".sanity_checked.release.gff3")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".final_sanity_check.log")
	shell:
		"sanity_check {input[0]} > {output[0]} 2> {log}"

rule gmc_generate_mikado_stats:
	input:
		rules.gmc_final_sanity_check.output[0]
	output:
		rules.gmc_final_sanity_check.output[0] + ".mikado_stats.txt",
		rules.gmc_final_sanity_check.output[0] + ".mikado_stats.tsv"
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="util stats"),
	shell:
		"{params.program_call} {input} --tab-stats {output[1]} > {output[0]} && " + \
		"parse_mikado_stats {output[0]} > {output[0]}.summary"

rule gmc_collect_biotype_conf_stats:
	input:
		rules.gmc_final_sanity_check.output[0]
	output:
		rules.gmc_final_sanity_check.output[0] + ".biotype_conf.tsv"
	run:
		import csv
		from collections import Counter
		biotypes, transcripts = dict(), dict()
		for row in csv.reader(open(input[0]), delimiter="\t"):
			if not row[0].startswith("#"):
				if row[2].lower() in {"gene", "mrna", "ncrna"}:
					attrib = dict(item.split("=") for item in row[8].strip().split(";"))
					if row[2].lower() == "gene":
						biotypes[attrib["ID"]] = attrib["biotype"]
					else:
						transcripts[attrib["ID"]] = (biotypes.get(attrib["Parent"], None), attrib["confidence"])
		with open(output[0], "w") as tbl_out:
			print("#1.transcript_id", "#2.biotype", "#3.confidence", sep="\t", file=tbl_out)
			for tid, bt_conf in sorted(transcripts.items()):
				print(tid, *bt_conf, sep="\t", file=tbl_out)
		counts = Counter(transcripts.values())
		with open(output[0].replace(".tsv", ".summary"), "w") as summary_out:
			for bt_conf, count in sorted(counts.items()):
				print(*bt_conf, count, sep="\t", file=summary_out)


rule gmc_extract_final_sequences:
	input:
		gff = rules.gmc_final_sanity_check.output[0],
		refseq = config["reference-sequence"]
	output:
		cdna = rules.gmc_final_sanity_check.output[0] + ".cdna.fasta",
		cds = rules.gmc_final_sanity_check.output[0] + ".cds.fasta",
		pep = rules.gmc_final_sanity_check.output[0] + ".pep.raw.fasta",
		tbl = rules.gmc_final_sanity_check.output[0] + ".gffread.table.txt"
	params:
		program_call = config["program_calls"]["gffread"],
		program_params = config["params"]["gffread"]["final"]
	shell:
		"{params.program_call} {input.gff} -g {input.refseq} {params.program_params} -w {output.cdna} -x {output.cds} -y {output.pep} -o {output.tbl}"

rule gmc_cleanup_final_proteins:
	input:
		rules.gmc_extract_final_sequences.output[2]
	output:
		rules.gmc_extract_final_sequences.output[2].replace(".raw.fasta", ".fasta")
	log:
		os.path.join(LOG_DIR, "cleanup_proteins.log")
	params:
		prefix = rules.gmc_extract_final_sequences.output[2].replace(".raw.fasta", ""),
		program_call = config["program_calls"]["prinseq"],
		program_params = config["params"]["prinseq"]
	shell:
		"{params.program_call} -fasta {input} {params.program_params} -out_good {params.prefix} -out_bad {params.prefix}.bad"

rule gmc_protein_completeness:
	input:
		proteins1 = rules.gmc_gffread_extract_proteins_post_pick.output[0],
		proteins2 = rules.gmc_extract_final_sequences.output[2]		
	output:
		tsv1 = os.path.join(config["outdir"], "mikado.annotation.protein_status.tsv"),
		summary1 = os.path.join(config["outdir"], "mikado.annotation.protein_status.summary"),
		tsv2 = rules.gmc_extract_final_sequences.output[2].replace(".fasta", ".protein_status.tsv"),
		summary2 = rules.gmc_extract_final_sequences.output[2].replace(".fasta", ".protein_status.summary")
	params:
		outdir = config["outdir"],
		prefix1 = "mikado.annotation",
		prefix2 = os.path.basename(rules.gmc_extract_final_sequences.output[2].replace(".fasta", ""))
	shell:
		"protein_completeness -o {params.outdir} -p {params.prefix1} {input.proteins1} && "  + \
		"protein_completeness -o {params.outdir} -p {params.prefix2} {input.proteins2}"
