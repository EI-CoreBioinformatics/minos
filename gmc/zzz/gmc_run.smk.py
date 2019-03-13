import os
import sys


EXTERNAL_METRICS_DIR = os.path.join(config["outdir"], "generate_metrics")
LOGDIR = os.path.join(config["outdir"], "logs")


def get_rnaseq(wc):
	# return " ".join(" ".join(pair) for pair in config["expression-runs"][wc.run])
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
		print("GENERATING TARGET:", run)
		OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "kallisto", run, "abundance.tsv"))

#if config.get("transcript-runs"):
#	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado", "vs_all", "MIKADO_DONE"))

for run in config.get("data", dict()).get("protein-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "proteins", run, "MIKADO_DONE"))
for run in config.get("data", dict()).get("transcript-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "transcripts", run, "MIKADO_DONE"))
for run in config.get("data", dict()).get("protein-seqs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], run, run + ".{}.tsv".format(config["blast-mode"])))

localrules: all, test_imitate_blast, gmc_metrics_blastp_combine

rule all:
	input:
		OUTPUTS
		# dynamic(os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt"))
		# dynamic(os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt.link"))
		#dynamic(os.path.join(config["outdir"], "tmp", "chunk-{chunk}.{run}.blastp.tsv"))
		

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
		os.path.join(LOGDIR, config["prefix"] + ".mikado_prepare.log")
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
		os.path.join(LOGDIR, config["prefix"] + ".gffread_extract.log")
	threads:
		1
	params:
		extract = "-y" if config["blast-mode"] == "blastp" else "-W -x"
	shell:
		"set +u && source cufflinks-2.2.1_gk && " + \
		"gffread {input.gtf} -g {input.refseq} {params.extract} {output[0]}.raw &> {log} && " + \
		#"sed 's/\.$/$/g' {output[0]}.raw > {output[0]}"
		"awk '/^[^>]/ {{ $1=gensub(\"\\\\.\", \"\", \"g\", $1) }} {{ print $0 }}' {output[0]}.raw > {output[0]}"






rule gmc_metrics_cpc2:
	input:
		rules.gmc_mikado_prepare.output[0]
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "CPC-2.0_beta", os.path.basename(rules.gmc_mikado_prepare.output[0]) + ".cpc2output.txt")
	params:
		x = 1
	log:
		os.path.join(LOGDIR, config["prefix"] + ".CPC2.log")
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
			os.path.join(LOGDIR, config["prefix"] + ".kallisto_index.log")
		shell:
			"set +u && source kallisto-0.44.0 && /usr/bin/time -v kallisto index -i {output[0]} {input[0]} &> {log}"

	rule gmc_metrics_kallisto_quant:
		input:
			index = rules.gmc_metrics_kallisto_index.output[0],
			reads = get_rnaseq
		output:
			os.path.join(EXTERNAL_METRICS_DIR, "kallisto", "{run}", "abundance.tsv")
		log:
			os.path.join(LOGDIR, config["prefix"] + ".{run}.kallisto.log")
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
		os.path.join(LOGDIR, config["prefix"] + ".mikado_compare.tran.{run}.log")
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
		os.path.join(LOGDIR, config["prefix"] + ".mikado_compare.prot.{run}.log")
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
		os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], "{run}", "{run}")
	log:
		os.path.join(LOGDIR, config["prefix"] + ".makeblastdb.{run}.log")
	params:
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], wildcards.run),
		db_prefix = lambda wildcards: wildcards.run
	shell:
		"set +u && source blast-2.6.0 && " + \
		"makeblastdb -in {input[0]} -out {params.outdir}/{params.db_prefix} -logfile {log} -dbtype prot && " + \
		"touch {output[0]}"


checkpoint gmc_chunk_proteins:
	input:
		# rules.gmc_gffread_generate_cds.output[0]
		rules.gmc_gffread_extract_proteins.output[0]
	output:
		# dynamic(os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt"))
		chunk_dir = directory(os.path.join(config["outdir"], "tmp"))
	log:
		os.path.join(LOGDIR, os.path.basename(rules.gmc_gffread_extract_proteins.output[0]) + ".chunk.log")
	params:
		chunksize = 1000,
		outdir = os.path.join(config["outdir"], "tmp")
	shell:
		"mkdir -p {params.outdir} && " + \
		# awk script by Pierre Lindenbaum https://www.biostars.org/p/13270/
		"awk 'BEGIN {{n=0;m=1;}} /^>/ {{ if (n%{params.chunksize}==0) {{f=sprintf(\"{params.outdir}/chunk-%d.txt\",m); m++;}}; n++; }} {{ print >> f }}' {input[0]} &> {log}"		

rule gmc_metrics_blastp_chunked:
	input:
		chunk = os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt"),
		db = rules.gmc_metrics_blastp_mkdb.output[0]
	output:
		os.path.join(config["outdir"], "tmp", "{run}", "chunk-{chunk}." + config["blast-mode"] + ".tsv")
	log:
		os.path.join(LOGDIR, "blast_logs", "chunk-{chunk}.{run}." + config["blast-mode"] + ".log")
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
	print("CHECKPOINT:", checkpoint_output)
	print(expand(
		os.path.join(config["outdir"], "tmp", "{run}", "chunk-{chunk}." + config["blast-mode"] + ".tsv"), #"post/{sample}/{i}.txt",
		run=wildcards.run,
		chunk=glob_wildcards(os.path.join(checkpoint_output, "chunk-{chunk}.txt")).chunk
	))
	return expand(
		os.path.join(config["outdir"], "tmp", "{run}", "chunk-{chunk}." + config["blast-mode"] + ".tsv"), #"post/{sample}/{i}.txt",
		run=wildcards.run,
		chunk=glob_wildcards(os.path.join(checkpoint_output, "chunk-{chunk}.txt")).chunk
	)


rule gmc_metrics_blastp_combine:
	input:
		aggregate_input # dynamic(os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt.blastp.tsv"))
	output:
		os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], "{run}", "{run}." + config["blast-mode"] + ".tsv")		
	shell:
		"cat {input} > {output[0]}"
