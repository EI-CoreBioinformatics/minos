import os
import sys


EXTERNAL_METRICS_DIR = os.path.join(config["outdir"], "generate_metrics")
LOG_DIR = os.path.join(config["outdir"], "logs")
TEMP_DIR = os.path.join(config["outdir"], "tmp")

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


for run in config.get("data", dict()).get("expression-runs", dict()):
	if config["use-tpm-for-picking"]:
		OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "kallisto", run, "abundance.tsv"))


for run in config.get("data", dict()).get("protein-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "proteins", run, "MIKADO_DONE"))
for run in config.get("data", dict()).get("transcript-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "transcripts", run, "MIKADO_DONE"))
for run in config.get("data", dict()).get("protein-seqs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], run, run + ".{}.tsv.tophit".format(config["blast-mode"])))

localrules: all, gmc_metrics_blastp_combine, gmc_metrics_generate_metrics_info

rule all:
	input:
		OUTPUTS,
		os.path.join(EXTERNAL_METRICS_DIR, "metrics_info.txt")
		
		

rule gmc_mikado_prepare:
	input:
		config["mikado-config-file"]
	output:
		os.path.join(config["outdir"], "mikado_prepared.fasta"),
		os.path.join(config["outdir"], "mikado_prepared.gtf")
	params:
		mikado = config["mikado-container"] + " mikado prepare",
		min_length = 100,
		outdir = config["outdir"]
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_prepare.log")
	threads:
		30
	shell:
		"singularity exec {params.mikado} --minimum_length {params.min_length} --json-conf {input[0]} --procs {threads} -od {params.outdir} &> {log}" 

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
		extract = "-y" if config["blast-mode"] == "blastp" else "-W -x"
	shell:
		"set +u && source cufflinks-2.2.1_gk && " + \
		"gffread {input.gtf} -g {input.refseq} {params.extract} {output[0]}.raw &> {log} && " + \
		"awk '/^[^>]/ {{ $1=gensub(\"\\\\.\", \"\", \"g\", $1) }} {{ print $0 }}' {output[0]}.raw > {output[0]}"


rule gmc_metrics_cpc2:
	input:
		rules.gmc_mikado_prepare.output[0]
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "CPC-2.0_beta", os.path.basename(rules.gmc_mikado_prepare.output[0]) + ".cpc2output.txt")
	params:
		x = 1
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".CPC2.log")
	threads:
		4
	shell:
		"set +u && source CPC-2.0_beta_py3_cs && " + \
		"/usr/bin/time -v CPC2.py -r -i {input[0]} -o {output[0]} &> {log} "


if config["use-tpm-for-picking"]:

	rule gmc_metrics_kallisto_index:
		input:
			rules.gmc_mikado_prepare.output[0]
		output:
			os.path.join(EXTERNAL_METRICS_DIR, "kallisto", os.path.basename(rules.gmc_mikado_prepare.output[0]) + ".idx")
		log:
			os.path.join(LOG_DIR, config["prefix"] + ".kallisto_index.log")
		shell:
			"set +u && source kallisto-0.44.0 && /usr/bin/time -v kallisto index -i {output[0]} {input[0]} &> {log}"

	rule gmc_metrics_kallisto_quant:
		input:
			index = rules.gmc_metrics_kallisto_index.output[0],
			reads = get_rnaseq
		output:
			os.path.join(EXTERNAL_METRICS_DIR, "kallisto", "{run}", "abundance.tsv")
		log:
			os.path.join(LOG_DIR, config["prefix"] + ".{run}.kallisto.log")
		params:
			stranded = lambda wildcards: "" if wildcards.run.endswith("_xx") else "--" + wildcards.run[-2:] + "-stranded",  
			bootstrap = 100,
			outdir = os.path.join(EXTERNAL_METRICS_DIR, "kallisto", "{run}")
		threads:
			32
		shell:
			"set +u && source kallisto-0.44.0 && /usr/bin/time -v kallisto quant {params.stranded} -i {input.index} -o {params.outdir} -b {params.bootstrap} --threads {threads} {input.reads} &> {log}"


rule gmc_metrics_mikado_compare_vs_transcripts:
	input:
		mika = rules.gmc_mikado_prepare.output[1],
		transcripts = get_transcript_alignments
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "transcripts", "{run}", "MIKADO_DONE")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_compare.tran.{run}.log")
	params:		                                                         	
		mikado = config["mikado-container"] + " mikado compare",
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "vs_transcripts", wildcards.run),
		transcripts = lambda wildcards: wildcards.run
	shell:
		"set +u && mkdir -p {params.outdir} && " + \
		"singularity exec {params.mikado} --extended-refmap -r {input.mika} -p {input.transcripts} -o {params.outdir}/vs_transcripts_{params.transcripts} &> {log} && " + \
		"touch {output[0]}"

		
rule gmc_metrics_mikado_compare_vs_proteins:
	input:
		mika = rules.gmc_mikado_prepare.output[1],
		proteins = get_protein_alignments
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "proteins", "{run}", "MIKADO_DONE")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_compare.prot.{run}.log")
	params:		                                                         	
		mikado = config["mikado-container"] + " mikado compare",
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "vs_proteins", wildcards.run),
		proteins = lambda wildcards: wildcards.run
	shell:
		"set +u && mkdir -p {params.outdir} && " + \
		"singularity exec {params.mikado} --exclude-utr --extended-refmap -r {input.mika} -p {input.proteins} -o {params.outdir}/vs_proteins_{params.proteins} &> {log} && " + \
		"touch {output[0]}"


rule gmc_metrics_blastp_mkdb:
	input:
		get_protein_sequences
	output:
		os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], "{run}", "blastdb", "{run}")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".makeblastdb.{run}.log")
	params:
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], wildcards.run, "blastdb"),
		db_prefix = lambda wildcards: wildcards.run
	shell:
		"set +u && source blast-2.6.0 && " + \
		"makeblastdb -in {input[0]} -out {params.outdir}/{params.db_prefix} -logfile {log} -dbtype prot && " + \
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
		blast = config["blast-mode"]
	shell:
		"set +u && source blast-2.6.0 && " + \
		"{params.blast} -query {input.chunk} -out {output[0]} -max_target_seqs 1 -evalue 1e-5 -num_threads {threads} " + \
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
		parse_script = "/ei/workarea/group-ga/Scripts/parse_blast_tabular_format_v0.4.pl",
		extract_script = "/ei/workarea/group-ga/Scripts/get_topHit.pl"
	shell:
		"set +u && source perl-5.20.1_gk && " + \
		"{params.parse_script} {input[0]} 0 0 query | {params.extract_script} - 1 > {output[0]}"


rule gmc_metrics_generate_metrics_info:
	input:
		OUTPUTS		
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "metrics_info.txt")
	shell:
		"ls {input} > {output}"
