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

# gmc_run4/generate_metrics/mikado_compare/vs_transcripts/Scallop_Old_leaf_transcripts/vs_transcripts_Scallop_Old_leaf_transcripts.refmap
for run in config.get("data", dict()).get("protein-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "proteins", run, "mikado_" + run + ".refmap"))
for run in config.get("data", dict()).get("transcript-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "transcripts", run, "mikado_" + run + ".refmap"))
for run in config.get("data", dict()).get("protein-seqs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], run, run + ".{}.tsv.tophit".format(config["blast-mode"])))


localrules: all, gmc_metrics_blastp_combine, gmc_metrics_generate_metrics_info, gmc_metrics_generate_metrics_matrix, gmc_parse_mikado_pick

rule all:
	input:
		OUTPUTS,
		os.path.join(EXTERNAL_METRICS_DIR, "metrics_info.txt"),
		os.path.join(config["outdir"], "MIKADO_SERIALISE_DONE"),
		os.path.join(config["outdir"], "mikado.subloci.gff3"),
		os.path.join(config["outdir"], "mikado.loci.gff3")
		

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

rule gmc_mikado_compare_index_reference:
	input:
		rules.gmc_mikado_prepare.output[1]
	output:
		rules.gmc_mikado_prepare.output[1] + ".midx"
	params:
		mikado = config["mikado-container"] + " mikado compare"
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_compare_index_reference.log")
	shell:
		"singularity exec {params.mikado} --index -r {input} &> {log}"

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
		mika = rules.gmc_mikado_compare_index_reference.output[0],
		# mika = rules.gmc_mikado_prepare.output[1],
		transcripts = get_transcript_alignments
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "transcripts", "{run}", "mikado_{run}.refmap")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_compare.tran.{run}.log")
	params:		                                                         	
		mikado = config["mikado-container"] + " mikado compare",
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "transcripts", wildcards.run),
		transcripts = lambda wildcards: wildcards.run
	shell:
		"set +u && mkdir -p {params.outdir} && " + \
		"singularity exec {params.mikado} --extended-refmap -r {input.mika} -p {input.transcripts} -o {params.outdir}/mikado_{params.transcripts} &> {log} && " + \
		"touch {output[0]}"

		
rule gmc_metrics_mikado_compare_vs_proteins:
	input:
		mika = rules.gmc_mikado_compare_index_reference.output[0],
		# mika = rules.gmc_mikado_prepare.output[1],
		proteins = get_protein_alignments
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "proteins", "{run}", "mikado_{run}.refmap")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_compare.prot.{run}.log")
	params:		                                                         	
		mikado = config["mikado-container"] + " mikado compare",
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "proteins", wildcards.run),
		proteins = lambda wildcards: wildcards.run
	shell:
		# --exclude-utr
		"set +u && mkdir -p {params.outdir} && " + \
		"singularity exec {params.mikado} --extended-refmap -r {input.mika} -p {input.proteins} -o {params.outdir}/mikado_{params.proteins} &> {log} && " + \
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
				if os.path.basename(cwd) == "CPC-2.0_beta":
					mclass, mid, path = "cpc", "cpc", os.path.join(cwd, files[0])
					rows.append((2, mclass, mid, path))
				elif os.path.basename(cwd) in {"proteins", "transcripts"}:
					mclass = "mikado"
					for mid in dirs:
						cwd, _, files = next(walk)
						path = glob.glob(os.path.join(cwd, "*.refmap"))[0]
						rows.append((0, mclass, mid, path))
				elif os.path.basename(cwd) in {"blastp", "blastx"}:
					mclass = "blast"
					for mid in dirs:
						cwd, _, files = next(walk)
						path = glob.glob(os.path.join(cwd, "*.tophit"))[0]
						rows.append((1, mclass, mid, path))
						# last block is not necessary if we clean up the blast databases before
						try:
							_ = next(walk)
						except StopIteration:
							pass

			for mid in config["data"].get("repeat-data", dict()):
				mclass = "repeat"
				path = config["data"]["repeat-data"][mid][0][0]
				rows.append((3, mclass, mid, path))

			for _, mclass, mid, path in sorted(rows, key=lambda x:x[0]):
			
				print(mclass, mid, path, sep="\t", file=sys.stderr)
				print(mclass, mid, os.path.abspath(path), sep="\t", file=metrics_info)


rule gmc_metrics_generate_metrics_matrix:
	input:
		rules.gmc_metrics_generate_metrics_info.output[0]
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "metrics_matrix.txt")
	params:
		#script = "/ei/workarea/group-ga/Scripts/create_scoring_metrics_for_Quillaja_saponaria_v0.1.pl"
		#script = "/ei/workarea/group-ga/Scripts/create_scoring_metrics_for_Melia_azedarach_v0.1.pl"
		#script = "/ei/workarea/group-pb/schudomc_sandbox/gmc_dev/testdata_MAZ/create_scoring_metrics.pl"
		script = "generate_metrics"
	log:
		os.path.join(LOG_DIR, "generate_metrics_matrix.log")
	shell:
		#"set +u && source perl-5.20.1_gk && " + \
		"{params.script} {input[0]} > {output[0]} 2> {log}"

"""
mikado serialise --procs 30 --json-conf Quillaja_saponaria.configuration.yaml --external-scores annotation_run1.metrics.txt
"""
rule gmc_mikado_serialise:
	input:
		config = config["mikado-config-file"],
		ext_scores = rules.gmc_metrics_generate_metrics_matrix.output[0],
		transcripts = rules.gmc_mikado_prepare.output[0]
	output:
		#os.path.join(config["outdir"], "mikado.subloci.gff3")
		os.path.join(config["outdir"], "MIKADO_SERIALISE_DONE")
	params:
		mikado = config["mikado-container"] + " mikado serialise",
		outdir = config["outdir"]
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_serialise.log")
	threads:
		32
	shell:
		"singularity exec {params.mikado} --transcripts {input.transcripts} --external-scores {input.ext_scores} --json-conf {input.config} --procs {threads} -od {params.outdir} &> {log} && " + \
		"touch {output[0]}"

rule gmc_mikado_pick:
	input:
		config = config["mikado-config-file"],
		gtf = rules.gmc_mikado_prepare.output[1],
		serialise_done = rules.gmc_mikado_serialise.output[0]
	output:
		loci = os.path.join(config["outdir"], "mikado.loci.gff3"),
		subloci = os.path.join(config["outdir"], "mikado.subloci.gff3")
	threads:
		30
	params:
		mikado = config["mikado-container"] + " mikado pick",
		outdir = config["outdir"]
	shell:
		"singularity exec {params.mikado} -od {params.outdir} --procs {threads} --json-conf {input.config} --subloci_out $(basename {output.subloci}) {input.gtf}"

rule gmc_parse_mikado_pick:
	input:
		loci = rules.gmc_mikado_pick.output[0]
	output:
		gff = os.path.join(config["outdir"], "Mikado.apollo.gff")
	shell:
		"parse_mikado_gff {input.loci} > {output.gff}"


#rule gmc_protein_completeness:
#	input:
			


