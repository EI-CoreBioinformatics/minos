import os
import sys


EXTERNAL_METRICS_DIR = os.path.join(config["outdir"], "generate_metrics")
LOGDIR = os.path.join(config["outdir"], "logs")


def get_rnaseq(wc):
	# return " ".join(" ".join(pair) for pair in config["expression-runs"][wc.run])
	return [item for sublist in config["expression-runs"][wc.run] for item in sublist]

def get_all_transcript_assemblies(wc):
	print([config["transcript-runs"][asm][0] for asm in config["transcript-runs"]])
	return [config["transcript-runs"][asm][0] for asm in config["transcript-runs"]]

def get_protein_alignments(wc):
	return config["protein-runs"][wc.run][0]
def get_transcript_alignments(wc):
	return config["transcript-runs"][wc.run][0]
def get_protein_sequences(wc):
	return config["protein-seqs"].get(wc.run, [""])[0]

OUTPUTS = [
	os.path.join(config["outdir"], "mikado_prepared.fasta"), 
	os.path.join(EXTERNAL_METRICS_DIR, "CPC-2.0_beta", "mikado_prepared.fasta.cpc2output.txt"),
	# os.path.join(config["outdir"], "logs", "mikado_prepared.fasta.leaff.log"),
	os.path.join(EXTERNAL_METRICS_DIR, "kallisto", "mikado_prepared.fasta.idx"),
	# os.path.join(LOGDIR, os.path.basename("mikado_prepared.gtf.cds.fasta.leaff.log"))
]

for run in config.get("expression-runs", {}):
	print("GENERATING TARGET:", run)
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "kallisto", run, "abundance.tsv"))

#if config.get("transcript-runs"):
#	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado", "vs_all", "MIKADO_DONE"))

for run in config.get("protein-runs", {}):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado", "proteins", run, "MIKADO_DONE"))
for run in config.get("transcript-runs", {}):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado", "transcripts", run, "MIKADO_DONE"))
for run in config.get("protein-seqs", {}):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "blastx", run, run + ".blastx.tsv"))

localrules: all, test_imitate_blast, gmc_metrics_blastx_combine

rule all:
	input:
		OUTPUTS
		# dynamic(os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt"))
		# dynamic(os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt.link"))
		#dynamic(os.path.join(config["outdir"], "tmp", "chunk-{chunk}.{run}.blastx.tsv"))
		

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

rule gmc_gffread_generate_cds:
	input:
		gtf = rules.gmc_mikado_prepare.output[1],
		refseq = config["reference-sequence"]
	output:
		rules.gmc_mikado_prepare.output[1] + ".cds.fasta"
	log:
		os.path.join(LOGDIR, config["prefix"] + ".gffread_cds.log")
	threads:
		1
	shell:
		"set +u && source cufflinks-2.2.1_gk && " + \
		"gffread {input.gtf} -g {input.refseq} -W -w {output[0]} &> {log}"

#rule gmc_leaff_chunk_cds:
#	input:
#		rules.gmc_gffread_generate_cds.output[0]
#	output:
#		dynamic(os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt"))
#	#log:
#	#	os.path.join(LOGDIR, os.path.basename(rules.gmc_gffread_generate_cds.output[0]) + ".leaff.log")
#	params:
#		chunksize = 1000,
#		outdir = os.path.join(config["outdir"], "tmp"),
#		wrapper = "/ei/workarea/group-ga/Scripts/leaff_v0.2.pl"
#	shell:
#		"set +u && source perl-5.20.1_gk && mkdir -p {params.outdir} && " + \
#        "{params.wrapper} {input[0]} {params.outdir}/chunk {params.chunksize} &> {output[0]}"


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

#rule gmc_metrics_leaff:
#	input:
#		rules.gmc_mikado_prepare.output[0]
#	output:
#		os.path.join(LOGDIR, os.path.basename(rules.gmc_mikado_prepare.output[0]) + ".leaff.log")
#	params:
#		chunksize = 500,
#		outdir = os.path.join(config["outdir"], "tmp"),
#		wrapper = "/ei/workarea/group-ga/Scripts/leaff_v0.2.pl"
#	shell:
#		"set +u && source perl-5.20.1_gk && mkdir -p {params.outdir} && " + \
#		"{params.wrapper} {input[0]} {params.outdir}/chunk {params.chunksize} &> {output[0]}"

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
		stranded = lambda wildcards: "" if wildcards.run.endswith("_xx") else "--" + wildcards.run[-2:] + "-stranded",  #if True else "", # TODO!
		bootstrap = 100,
		outdir = os.path.join(EXTERNAL_METRICS_DIR, "kallisto", "{run}")
	threads:
		32
	shell:
		"set +u && source kallisto-0.44.0 && /usr/bin/time -v kallisto quant {params.stranded} -i {input.index} -o {params.outdir} -b {params.bootstrap} --threads {threads} {input.reads} &> {log}"

"""
rule gmc_metrics_mikado_vs_all:
	input:
		mika = rules.gmc_mikado_prepare.output[1],
		tr_assemblies = get_all_transcript_assemblies
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "mikado", "vs_all", "MIKADO_DONE")
	log:
		os.path.join(LOGDIR, config["prefix"] + ".mikado.vs_all.log")
	params:		
		mikado = config["mikado-container"] + " mikado compare",
		outdir = os.path.join(EXTERNAL_METRICS_DIR, "mikado", "vs_all")		
	shell:
		"set +u && mkdir -p {params.outdir} && cat {input.tr_assemblies} > {params.outdir}/gmc_all_transcripts_merged.gtf && " + \
		"singularity exec {params.mikado} --extended-refmap -r {input.mika} -p {params.outdir}/gmc_all_transcripts_merged.gtf -o {params.outdir}/vs_all &> {log} && " + \
		"touch {output[0]}"
"""

rule gmc_metrics_mikado_vs_transcripts:
	input:
		mika = rules.gmc_mikado_prepare.output[1],
		transcripts = get_transcript_alignments
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "mikado", "transcripts", "{run}", "MIKADO_DONE")
	log:
		os.path.join(LOGDIR, config["prefix"] + ".mikado_compare.tran.{run}.log")
	params:		                                                         	
		mikado = config["mikado-container"] + " mikado compare",
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, "mikado", "vs_transcripts", wildcards.run),
		transcripts = lambda wildcards: wildcards.run
	shell:
		"set +u && mkdir -p {params.outdir} && " + \
		"singularity exec {params.mikado} --extended-refmap -r {input.mika} -p {input.transcripts} -o {params.outdir}/vs_transcripts_{params.transcripts} &> {log} && " + \
		"touch {output[0]}"

		
rule gmc_metrics_mikado_vs_proteins:
	input:
		mika = rules.gmc_mikado_prepare.output[1],
		proteins = get_protein_alignments
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "mikado", "proteins", "{run}", "MIKADO_DONE")
	log:
		os.path.join(LOGDIR, config["prefix"] + ".mikado_compare.prot.{run}.log")
	params:		                                                         	
		mikado = config["mikado-container"] + " mikado compare",
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, "mikado", "vs_proteins", wildcards.run),
		proteins = lambda wildcards: wildcards.run
	shell:
		"set +u && mkdir -p {params.outdir} && " + \
		"singularity exec {params.mikado} --exclude-utr --extended-refmap -r {input.mika} -p {input.proteins} -o {params.outdir}/vs_proteins_{params.proteins} &> {log} && " + \
		"touch {output[0]}"

rule gmc_metrics_blastx_mkdb:
	input:
		get_protein_sequences
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "blastx", "{run}", "{run}")
	log:
		os.path.join(LOGDIR, config["prefix"] + ".makeblastdb.{run}.log")
	params:
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, "blastx", wildcards.run),
		db_prefix = lambda wildcards: wildcards.run
	shell:
		"set +u && source blast-2.6.0 && " + \
		"makeblastdb -in {input[0]} -out {params.outdir}/{params.db_prefix} -logfile {log} -dbtype prot && " + \
		"touch {output[0]}"

# rule test_imitate_blast:
# 	input:
# 		os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt")
# 	output:
# 		os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt.link")
# 	shell:
# 		"ln -sf $(basename {input[0]}) {output[0]}"

checkpoint gmc_leaff_chunk_cds:
	input:
		rules.gmc_gffread_generate_cds.output[0]
	output:
		# dynamic(os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt"))
		chunk_dir = directory(os.path.join(config["outdir"], "tmp"))
	log:
		os.path.join(LOGDIR, os.path.basename(rules.gmc_gffread_generate_cds.output[0]) + ".leaff.log")
	params:
		chunksize = 1000,
		outdir = os.path.join(config["outdir"], "tmp"),
		wrapper = "/ei/workarea/group-ga/Scripts/leaff_v0.2.pl"
	shell:
		"set +u && source perl-5.20.1_gk && mkdir -p {params.outdir} && " + \
        "{params.wrapper} {input[0]} {params.outdir}/chunk {params.chunksize} &> {log}"

rule gmc_metrics_blastx_chunked:
	input:
		chunk = os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt"),
		db = rules.gmc_metrics_blastx_mkdb.output[0]
	output:
		os.path.join(config["outdir"], "tmp", "{run}", "chunk-{chunk}.blastx.tsv")
	log:
		os.path.join(LOGDIR, "blastx_logs", "chunk-{chunk}.{run}.blastx.log")
	threads:
		8
	shell:
		"set +u && source blast-2.6.0 && " + \
		"blastx -query {input.chunk} -out {output[0]} -max_target_seqs 1 -evalue 1e-5 -num_threads {threads} " + \
		"-db {input.db} -outfmt \"6 qseqid sseqid pident qstart qend sstart send qlen slen length nident mismatch positive gapopen gaps evalue bitscore\" &> {log}"

def aggregate_input(wildcards):
	checkpoint_output = checkpoints.gmc_leaff_chunk_cds.get(**wildcards).output.chunk_dir
	print("CHECKPOINT:", checkpoint_output)
	print(expand(
		os.path.join(config["outdir"], "tmp", "{run}", "chunk-{chunk}.blastx.tsv"), #"post/{sample}/{i}.txt",
		run=wildcards.run,
		chunk=glob_wildcards(os.path.join(checkpoint_output, "chunk-{chunk}.txt")).chunk
	))
	return expand(
		os.path.join(config["outdir"], "tmp", "{run}", "chunk-{chunk}.blastx.tsv"), #"post/{sample}/{i}.txt",
		run=wildcards.run,
		chunk=glob_wildcards(os.path.join(checkpoint_output, "chunk-{chunk}.txt")).chunk
	)


rule gmc_metrics_blastx_combine:
	input:
		aggregate_input # dynamic(os.path.join(config["outdir"], "tmp", "chunk-{chunk}.txt.blastx.tsv"))
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "blastx", "{run}", "{run}.blastx.tsv")		
	shell:
		"cat {input} > {output[0]}"
