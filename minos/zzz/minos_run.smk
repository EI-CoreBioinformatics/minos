import os
import sys
import pathlib
from minos.minos_configure import ExternalMetrics
from eicore.hpc_config import HpcConfig
HPC_CONFIG = HpcConfig(config["hpc_config"])

EXTERNAL_METRICS_DIR = os.path.join(config["outdir"], "generate_metrics")
pathlib.Path(EXTERNAL_METRICS_DIR).mkdir(exist_ok=True, parents=True)

LOG_DIR = os.path.join(config["outdir"], "logs")
TEMP_DIR = os.path.join(config["outdir"], "tmp")
AUGUSTUS_CONFIG_DATA = config["paths"]["augustus_config_data"]

# check if we have RNA-Seq expression data
EXPRESSION_OK = len(config["data"]["expression-runs"])

def get_rnaseq(wc):
	return [item for sublist in config["data"]["expression-runs"][wc.run] for item in sublist]

def get_all_transcript_assemblies(wc):
	return [config["data"]["transcript-runs"][asm][0] for asm in config["data"]["transcript-runs"]]

def get_protein_alignments(wc):
	return config["data"]["protein-runs"][wc.run][0]
def get_transcript_alignments(wc):
	return config["data"]["transcript-runs"][wc.run][0]
def get_protein_sequences(wc):
	return config["data"]["protein-seqs"].get(wc.run, [""])[0]
def get_repeat_data(wc):
	return config["data"]["repeat-data"][wc.run][0]
def get_transcript_models(wc):
	return config["data"]["transcript_models"][wc.run]

# output management
POST_PICK_PREFIX = "mikado.annotation"
RELEASE_PREFIX = config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1")
RESULTS_DIR = os.path.join(config["outdir"], "results")

OUTPUTS = [
	os.path.join(config["outdir"], "mikado_prepared.fasta"),
	os.path.join(config["outdir"], "mikado_prepared.exon.gff"),
	os.path.join(EXTERNAL_METRICS_DIR, "CPC-2.0_beta", "mikado_prepared.fasta.cpc2output.txt"),
]

if EXPRESSION_OK > 0:
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "kallisto", "mikado_prepared.fasta.idx"))

POST_PICK_EXPRESSION = list()
for run in config.get("data", dict()).get("expression-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "kallisto", run, "abundance.tsv"))
	POST_PICK_EXPRESSION.append(os.path.join(config["outdir"], "kallisto", run, "abundance.tsv"))

for run in config.get("data", dict()).get("protein-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "proteins", run, "mikado_" + run + ".refmap"))
for run in config.get("data", dict()).get("transcript-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "transcripts", run, "mikado_" + run + ".refmap"))
for run in config.get("data", dict()).get("protein-seqs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], run, run + ".DIAMOND.{}.tsv.tophit".format(config["blast-mode"]) if config["use-diamond"] else run + ".BLAST.{}.tsv.tophit".format(config["blast-mode"])))
for run in config.get("data", dict()).get("repeat-data", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "repeats", "{}.no_strand.exon.gff".format(run)))
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "repeats", "{}.no_strand.exon.gff.cbed".format(run)))
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "repeats", "{}.no_strand.exon.gff.cbed.parsed.txt".format(run)))


BUSCO_PRECOMPUTED = """
	cp {input[0]} {output[0]} && cp {input[1]} {output[1]} && cp {input[2]} {output[2]}
""".strip().replace("\n\t", " ")

BUSCO_CMD = """
	mkdir -p {BUSCO_PATH}/logs
	&& cfgdir={BUSCO_PATH}/config/{params.busco_stage}/{params.run} && mkdir -p $cfgdir
	&& {params.copy_busco} {AUGUSTUS_CONFIG_DATA} $cfgdir/
	&& export AUGUSTUS_CONFIG_PATH=$cfgdir/config
	&& rm -rvf {BUSCO_PATH}/runs/{params.busco_stage}/{params.run}
	&& mkdir -p {BUSCO_PATH}/runs/{params.busco_stage}/{params.run}
	&& cd {BUSCO_PATH}/runs/{params.busco_stage}/{params.run}
	&& {params.program_call} {params.program_params} -i {params.input} -c {threads} -m {params.busco_mode} --force -l {params.lineage_path} -o {params.run} &> {log}
	&& mv {params.run}/* . && rm -rf {params.run}
	&& rm -rf run_*/hmmer_output
	&& touch {output[2]}
	&& rm -rf $cfgdir
""".strip().replace("\n\t", " ")

TX2GENE_MAPS = [
	os.path.join(config["outdir"], "tx2gene", tm + ".tx2gene")
	for tm in config["data"]["transcript_models"]
]

BUSCO_PATH = os.path.abspath(os.path.join(config["outdir"], "busco"))
BUSCO_LINEAGE = os.path.basename(config["busco_analyses"]["lineage"]) if config["busco_analyses"]["lineage"] is not None else None

BUSCO_ANALYSES = list()
BUSCO_PROTEIN_PREPARE_RUNS = list()
if config["busco_analyses"]["proteins"]:
	BUSCO_PROTEIN_PREPARE_RUNS.extend(
		os.path.join(BUSCO_PATH, "runs", "proteins_prepare", "{tm}", "run_{lineage}", busco_file).format(tm=tm, lineage=BUSCO_LINEAGE)
		for tm in config["data"]["transcript_models"] for busco_file in ("short_summary.txt", "full_table.tsv", "missing_busco_list.tsv")
	)
	BUSCO_ANALYSES.extend(
		os.path.join(BUSCO_PATH, "runs", "proteins_final", "proteins_final", "run_{lineage}", busco_file).format(lineage=BUSCO_LINEAGE)
		for busco_file in ("short_summary.txt", "full_table.tsv", "missing_busco_list.tsv")
	)

if config["busco_analyses"]["transcriptome"]:
	BUSCO_ANALYSES.extend(
		os.path.join(BUSCO_PATH, "runs", "transcripts_prepare", "{tm}", "run_{lineage}", busco_file).format(tm=tm, lineage=BUSCO_LINEAGE)
		for tm in config["data"]["transcript_models"] for busco_file in ("short_summary.txt", "full_table.tsv", "missing_busco_list.tsv")
	)
	BUSCO_ANALYSES.extend(
		os.path.join(BUSCO_PATH, "runs", "transcripts_final", "transcripts_final", "run_{lineage}", busco_file).format(lineage=BUSCO_LINEAGE)
		for busco_file in ("short_summary.txt", "full_table.tsv", "missing_busco_list.tsv")
	)

if config["busco_analyses"]["genome"] or config["busco_analyses"]["precomputed_genome"]:
	runid = "precomputed" if config["busco_analyses"]["precomputed_genome"] else BUSCO_LINEAGE
	BUSCO_ANALYSES.extend(
		os.path.join(BUSCO_PATH, "runs", "genome", "genome", "run_{}".format(runid), busco_file)
		for busco_file in ("short_summary.txt", "full_table.tsv", "missing_busco_list.tsv")
	)
	BUSCO_GENOME_OUTPUT_SUMMARY, BUSCO_GENOME_OUTPUT_FULL_TABLE, BUSCO_GENOME_MISSING = BUSCO_ANALYSES[-3:]


BUSCO_COPY = [
	(path, os.path.join(RESULTS_DIR, "busco", path.replace(os.path.join(BUSCO_PATH, "runs", ""), "")))
	for path in BUSCO_PROTEIN_PREPARE_RUNS + BUSCO_ANALYSES
]

BUSCO_COPY_SOURCES = [src for src, dest in BUSCO_COPY]
BUSCO_COPY_TARGETS = [dest for src, dest in BUSCO_COPY]

BUSCO_TABLE = ""
if BUSCO_ANALYSES or BUSCO_PROTEIN_PREPARE_RUNS:
	BUSCO_ANALYSES.extend(TX2GENE_MAPS)
	BUSCO_TABLE = os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".release.busco_final_table.tsv")
	OUTPUTS.extend(BUSCO_PROTEIN_PREPARE_RUNS)

CONFIG_FILES = [
	os.path.join(config["outdir"], "{prefix}." + config_file).format(prefix=prefix) for prefix in [config["prefix"]] for config_file in ("scoring.yaml", "run_config.yaml", "mikado_config.yaml")
]

CONFIG_COPY = [
	(path, os.path.join(RESULTS_DIR, "config", path.replace(os.path.join(config["outdir"], ""), "")))
	for path in CONFIG_FILES
]

CONFIG_COPY_SOURCES = [src for src, dest in CONFIG_COPY]
CONFIG_COPY_TARGETS = [dest for src, dest in CONFIG_COPY]

localrules:
	all,
	minos_extract_exons,
	minos_metrics_repeats_convert,
	minos_mikado_prepare_extract_coords,
	minos_mikado_pick_extract_coords,
	minos_metrics_parse_repeat_coverage,
	minos_metrics_blastp_combine,
	minos_metrics_generate_metrics_info,
	minos_parse_mikado_pick,
	minos_gffread_extract_sequences_post_pick,
	minos_gffread_extract_sequences,
	minos_gff_genometools_check_post_pick,
	minos_calculate_cds_lengths_post_pick,
	minos_extract_final_sequences,
	split_proteins_prepare,
	split_transcripts_prepare,
	busco_copy_results,
	config_copy_results,
	busco_concat_protein_metrics,
	busco_summary,
	minos_create_release_metrics,
	minos_collate_metric_oddities,
	minos_summarise_collapsed_metrics


rule all:
	input:
		OUTPUTS,
		os.path.join(EXTERNAL_METRICS_DIR, "metrics_info.txt"),
		os.path.join(config["outdir"], "MIKADO_SERIALISE_DONE"),
		os.path.join(config["outdir"], "mikado.subloci.gff3"),
		os.path.join(config["outdir"], "mikado.monoloci.gff3"),
		os.path.join(config["outdir"], "mikado.loci.gff3"),
		os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".release.metrics_oddities.tsv"),
		#os.path.join(RESULTS_DIR, "mikado.monoloci.metrics_oddities.tsv"),
		#os.path.join(RESULTS_DIR, "mikado.loci.metrics_oddities.tsv"),
		[
			os.path.join(config["outdir"], POST_PICK_PREFIX + suffix)
			for suffix in {".gff", ".proteins.fasta", ".table.txt", ".cds.fasta", ".cds.fasta.lengths", ".cdna.fasta"}
		],
		os.path.join(config["outdir"], "kallisto", POST_PICK_PREFIX + ".cdna.fasta.idx") if EXPRESSION_OK > 0 else [],
		POST_PICK_EXPRESSION,
		[
			os.path.join(config["outdir"], POST_PICK_PREFIX + suffix)
			for suffix in {".gt_checked.gff", ".collapsed_metrics.tsv", ".gt_checked.validation_report.txt", ".release.unsorted.gff3", ".release_browser.unsorted.gff3"}
		],
		os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".release.collapsed_metrics.summary.tsv"),
		[
			os.path.join(RESULTS_DIR, RELEASE_PREFIX + suffix)
			for suffix in {".release.gff3", ".release.browser.gff3",
							".sanity_checked.release.gff3", ".release.gff3.mikado_stats.txt", ".release.gff3.mikado_stats.tsv"}
		],
		[
			os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".release.gff3") + ".{}.fasta".format(dtype) for dtype in {"cdna", "cds", "pep.raw", "pep"}
		],
		[
			os.path.join(RESULTS_DIR, RELEASE_PREFIX + suffix)
			for suffix in {".release.gff3.final_table.tsv", ".release.gff3.biotype_conf.summary", ".release.metrics.tsv"}
		],
		BUSCO_ANALYSES + ([BUSCO_TABLE] if BUSCO_TABLE else []),
		BUSCO_COPY_TARGETS,
		CONFIG_COPY_TARGETS



rule minos_mikado_prepare:
	input:
		config["mikado-config-file"]
	output:
		os.path.join(config["outdir"], "mikado_prepared.fasta"),
		os.path.join(config["outdir"], "mikado_prepared.gtf"),
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="prepare"),
		program_params = config["params"]["mikado"]["prepare"],
		outdir = config["outdir"]
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_prepare.log")
	threads:
		HPC_CONFIG.get_cores("minos_mikado_prepare")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_mikado_prepare") * attempt
	shell:
		"{params.program_call} {params.program_params} --json-conf {input[0]} --procs {threads} -od {params.outdir} &> {log}"

rule minos_mikado_prepare_extract_coords:
	input:
		rules.minos_mikado_prepare.output[1]
	output:
		rules.minos_mikado_prepare.output[1] + ".coords"
	run:
		from minos.scripts.extract_coords import extract_coords
		extract_coords(input[0], output[0], filetype="gtf")

rule minos_generate_tx2gene_maps:
	input:
		get_transcript_models
	output:
		os.path.join(config["outdir"], "tx2gene", os.path.basename("{run}") + ".tx2gene")
	threads:
		HPC_CONFIG.get_cores("minos_generate_tx2gene_maps")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_generate_tx2gene_maps")
	run:
		from minos.scripts.generate_tx2gene_maps import Gxf8Parser
		with open(output[0], "w") as _out:
			for tid, gid in Gxf8Parser(input[0]).parse_file(input[0], wildcards.run):
				print(tid, gid, sep="\t", flush=True, file=_out)

rule minos_extract_exons:
	input:
		rules.minos_mikado_prepare.output[1]
	output:
		rules.minos_mikado_prepare.output[1].replace(".gtf", ".exon.gff")
	run:
		from minos.scripts.extract_exons import extract_exons
		extract_exons(input[0], output[0])

rule minos_mikado_compare_index_reference:
	input:
		rules.minos_mikado_prepare.output[1]
	output:
		rules.minos_mikado_prepare.output[1] + ".midx"
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="compare"),
		program_params = config["params"]["mikado"]["compare"]["index"]
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".mikado_compare_index_reference.log")
	threads:
		HPC_CONFIG.get_cores("minos_mikado_compare_index_reference")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_mikado_compare_index_reference") * attempt
	shell:
		"{params.program_call} {params.program_params} -r {input} &> {log}"

rule minos_gffread_extract_sequences:
	input:
		gtf = rules.minos_mikado_prepare.output[1],
		refseq = config["reference-sequence"]
	output:
		rules.minos_mikado_prepare.output[1] + (".prot.fasta" if config["blast-mode"] == "blastp" else ".cds.fasta"),
		rules.minos_mikado_prepare.output[1] + ".cdna.fasta"
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".gffread_extract.log")
	params:
		program_call = config["program_calls"]["gffread"],
		program_params = config["params"]["gffread"][config["blast-mode"]],
		output_params = "-W -x" if config["blast-mode"] == "blastx" else ("-y" if config["blast-mode"] == "blastp" else "")
	shell:
		"{params.program_call} {input.gtf} -g {input.refseq} {params.program_params} -W -w {output[1]} {params.output_params} {output[0]} &> {log}"

rule minos_metrics_repeats_convert:
	input:
		get_repeat_data
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "repeats", "{run}.converted.gff"),
		os.path.join(EXTERNAL_METRICS_DIR, "repeats", "{run}.no_strand.exon.gff")
	run:
		from minos.scripts.parse_repeatmasker import parse_repeatmasker
		parse_repeatmasker(input[0], output[0], output[1], wildcards.run)

rule minos_metrics_bedtools_repeat_coverage:
	input:
		rules.minos_metrics_repeats_convert.output[1],
		rules.minos_extract_exons.output[0]
	output:
		rules.minos_metrics_repeats_convert.output[1] + ".cbed",
	params:
		program_call = config["program_calls"]["bedtools"]["coverageBed"]
	threads:
		HPC_CONFIG.get_cores("minos_metrics_bedtools_repeat_coverage")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_metrics_bedtools_repeat_coverage") * attempt
	shell:
		"""
		{params.program_call} -a {input[1]} -b {input[0]} > {output[0]}.tmp
		&& cut -f 9 {output[0]}.tmp | cut -d ';' -f 1 | paste - {output[0]}.tmp | sort -k1,1V | cut -f 2- > {output[0]}
		&& rm {output[0]}.tmp
		""".strip().replace("\n\t", " ")


rule minos_metrics_parse_repeat_coverage:
	input:
		rules.minos_metrics_bedtools_repeat_coverage.output[0],
	output:
		rules.minos_metrics_bedtools_repeat_coverage.output[0] + ".parsed.txt",
	run:
		from minos.scripts.parse_cbed_stats import parse_cbed
		with open(output[0], "w") as outstream, open(input[0]) as instream:
			parse_cbed(instream, outstream=outstream)


rule minos_metrics_cpc2:
	input:
		rules.minos_mikado_prepare.output[0]
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "CPC-2.0_beta", os.path.basename(rules.minos_mikado_prepare.output[0]) + ".cpc2output.txt")
	params:
		program_call = config["program_calls"]["cpc2"],
		program_params = config["params"]["cpc2"]
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".CPC2.log")
	threads:
		HPC_CONFIG.get_cores("minos_metrics_cpc2")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_metrics_cpc2") * attempt
	shell:
		"{params.program_call} {params.program_params} -i {input[0]} -o {output[0]} &> {log}"


rule minos_metrics_kallisto_index:
	input:
		rules.minos_mikado_prepare.output[0]
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "kallisto", os.path.basename(rules.minos_mikado_prepare.output[0]) + ".idx")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".kallisto_index.log")
	params:
		program_call = config["program_calls"]["kallisto"].format(program="index"),
		program_params = config["params"].get("kallisto", {}).get("index", "")
	threads:
		HPC_CONFIG.get_cores("minos_metrics_kallisto_index")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_metrics_kallisto_index") * attempt
	shell:
		"{params.program_call} {params.program_params} -i {output[0]} {input[0]} &> {log}"

rule minos_metrics_kallisto_quant:
	input:
		index = rules.minos_metrics_kallisto_index.output[0],
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
		HPC_CONFIG.get_cores("minos_metrics_kallisto_quant")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_metrics_kallisto_quant") * attempt
	shell:
		"{params.program_call} {params.program_params} {params.stranded} -i {input.index} -o {params.outdir} --threads {threads} {input.reads} &> {log}"

rule minos_metrics_mikado_compare_vs_transcripts:
	input:
		midx = rules.minos_mikado_compare_index_reference.output[0],
		mika = rules.minos_mikado_prepare.output[1],
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
	threads:
		HPC_CONFIG.get_cores("minos_metrics_mikado_compare_vs_transcripts")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_metrics_mikado_compare_vs_transcripts") * attempt
	shell:
		"mkdir -p {params.outdir}" + \
		" && {params.program_call} {params.program_params} -r {input.mika} -p {input.transcripts} -o {params.outdir}/mikado_{params.transcripts} &> {log}" + \
		" && touch {output[0]}"

rule minos_metrics_mikado_compare_vs_proteins:
	input:
		midx = rules.minos_mikado_compare_index_reference.output[0],
		mika = rules.minos_mikado_prepare.output[1],
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
	threads:
		HPC_CONFIG.get_cores("minos_metrics_mikado_compare_vs_proteins")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_metrics_mikado_compare_vs_proteins") * attempt
	shell:
		"mkdir -p {params.outdir}" + \
		" && {params.program_call} {params.program_params} -r {input.mika} -p {input.proteins} -o {params.outdir}/mikado_{params.proteins} &> {log}" + \
		" && touch {output[0]}"

BLASTDB_CMD = "{params.program_call} {params.program_params} -in {input[0]} -out {params.outdir}/{params.db_prefix} -logfile {log}"
DIAMONDDB_CMD = "{params.program_call} {params.program_params} --in {input[0]} --db {params.outdir}/{params.db_prefix} --verbose &> {log}"

rule minos_metrics_blastp_mkdb:
	input:
		get_protein_sequences
	output:
		os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], "{run}", "diamonddb" if config["use-diamond"] else "blastdb", "{run}")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".makedb.{run}.log" if config["use-diamond"] else config["prefix"] + ".makeblastdb.{run}.log")
	params:
		program_call = config["program_calls"]["diamond"].format(program="makedb") if config["use-diamond"] else config["program_calls"]["blast"].format(program="makeblastdb"),
		program_params = config["params"]["diamond"]["makedb"] if config["use-diamond"] else config["params"]["blast"]["makeblastdb"],
		outdir = lambda wildcards: os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], wildcards.run, "diamonddb" if config["use-diamond"] else "blastdb"),
		db_prefix = lambda wildcards: wildcards.run
	threads:
		HPC_CONFIG.get_cores("minos_metrics_blastp_mkdb")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_metrics_blastp_mkdb") * attempt
	shell:
		"{cmd} && touch {{output[0]}}".format(cmd=DIAMONDDB_CMD if config["use-diamond"] else BLASTDB_CMD)

checkpoint minos_chunk_proteins:
	input:
		rules.minos_gffread_extract_sequences.output[0]
	output:
		chunk_dir = directory(os.path.join(TEMP_DIR, "chunked_proteins"))
	log:
		os.path.join(LOG_DIR, os.path.basename(rules.minos_gffread_extract_sequences.output[0]) + ".chunk.log")
	params:
		chunksize = 1000,
		outdir = os.path.join(TEMP_DIR, "chunked_proteins")
	threads:
		HPC_CONFIG.get_cores("minos_chunk_proteins")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_chunk_proteins") * attempt
	shell:
		"mkdir -p {params.outdir}" + \
		# awk script by Pierre Lindenbaum https://www.biostars.org/p/13270/
		" && awk 'BEGIN {{n=0;m=1;}} /^>/ {{ if (n%{params.chunksize}==0) {{f=sprintf(\"{params.outdir}/chunk-%d.txt\",m); m++;}}; n++; }} {{ print >> f }}' {input[0]} &> {log}"

BLAST_CMD = "{params.program_call} {params.program_params} -query {input.chunk} -out {output[0]} -num_threads {threads} " + \
			"-db {input.db} -outfmt \"6 qseqid sseqid pident qstart qend sstart send qlen slen length nident mismatch positive gapopen gaps evalue bitscore\" &> {log}"
DIAMOND_CMD = "{params.program_call} {params.program_params} --query {input.chunk} --out {output[0]} --threads {threads} " + \
			"--db {input.db} --outfmt 6 qseqid sseqid pident qstart qend sstart send qlen slen length nident mismatch positive gapopen gaps evalue bitscore &> {log}"

rule minos_metrics_blastp_chunked:
	input:
		chunk = os.path.join(TEMP_DIR, "chunked_proteins", "chunk-{chunk}.txt"),
		db = rules.minos_metrics_blastp_mkdb.output[0]
	output:
		os.path.join(TEMP_DIR, "chunked_proteins", "{run}", "chunk-{chunk}.DIAMOND." + config["blast-mode"] + ".tsv" if config["use-diamond"] else "chunk-{chunk}.BLAST." + config["blast-mode"] + ".tsv")
	log:
		os.path.join(LOG_DIR, "blast_logs", "chunk-{chunk}.{run}.DIAMOND." + config["blast-mode"] + ".log" if config["use-diamond"] else "chunk-{chunk}.{run}.BLAST." + config["blast-mode"] + ".log")
	params:
		program_call = config["program_calls"]["diamond" if config["use-diamond"] else "blast"].format(program=config["blast-mode"]),
		program_params = config["params"]["diamond" if config["use-diamond"] else "blast"][config["blast-mode"]]
	threads:
		HPC_CONFIG.get_cores("minos_metrics_blastp_chunked")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_metrics_blastp_chunked") * attempt
	shell:
		"{cmd}".format(cmd=DIAMOND_CMD if config["use-diamond"] else BLAST_CMD)

def aggregate_blastp_input(wildcards):
	checkpoint_output = checkpoints.minos_chunk_proteins.get(**wildcards).output.chunk_dir
	return expand(
		os.path.join(TEMP_DIR, "chunked_proteins", "{run}", "chunk-{chunk}.DIAMOND." + config["blast-mode"] + ".tsv" if config["use-diamond"] else "chunk-{chunk}.BLAST." + config["blast-mode"] + ".tsv"),
		run=wildcards.run,
		chunk=glob_wildcards(os.path.join(checkpoint_output, "chunk-{chunk}.txt")).chunk
	)


rule minos_metrics_blastp_combine:
	input:
		aggregate_blastp_input
	output:
		os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], "{run}", "{run}.DIAMOND." + config["blast-mode"] + ".tsv" if config["use-diamond"] else "{run}.BLAST." + config["blast-mode"] + ".tsv")
	threads:
		HPC_CONFIG.get_cores("minos_metrics_blastp_combine")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_metrics_blastp_combine") * attempt
	run:
		with open(output[0], "w") as blast_out:
			for f in input:
				print(open(f).read(), end="", flush=True, file=blast_out)
				os.remove(f)


rule minos_metrics_blastp_tophit:
	input:
		rules.minos_metrics_blastp_combine.output[0]
	output:
		rules.minos_metrics_blastp_combine.output[0] + ".tophit"
	params:
		pident_threshold = config["params"]["diamond" if config["use-diamond"] else "blast"]["tophit"]["pident_threshold"],
		qcov_threshold = config["params"]["diamond" if config["use-diamond"] else "blast"]["tophit"]["qcov_threshold"]
	threads:
		HPC_CONFIG.get_cores("minos_metrics_blastp_tophit")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_metrics_blastp_tophit") * attempt
	run:
		from minos.scripts.get_blast_tophit import get_blast_tophit
		get_blast_tophit(input[0], output[0], params.pident_threshold, params.qcov_threshold)

rule busco_concat_protein_metrics:
	input:
		BUSCO_PROTEIN_PREPARE_RUNS
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "busco_proteins", "busco_proteins.tsv")
	run:
		with open(output[0], "w") as concat_out:
			for f in input:
				if os.path.basename(f) == "full_table.tsv":
					print("#" + f, flush=True, file=concat_out)
					print(open(f).read(), end="", flush=True, file=concat_out)


rule minos_metrics_generate_metrics_info:
	input:
		prev_outputs = OUTPUTS,
		bproteins = os.path.join(EXTERNAL_METRICS_DIR, "busco_proteins", "busco_proteins.tsv")
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "metrics_info.txt")
	run:
		from minos.scripts.generate_metrics_info import generate_metrics_info
		busco_data = list(config["data"].get("busco-data", {"busco_proteins": ""}).keys())[0]
		generate_metrics_info(EXTERNAL_METRICS_DIR, output[0], busco_data, config["use-diamond"])


rule minos_metrics_generate_metrics_matrix:
	input:
		rules.minos_metrics_generate_metrics_info.output[0]
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "metrics_matrix.txt")
	log:
		os.path.join(LOG_DIR, "generate_metrics_matrix.log")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_metrics_generate_metrics_matrix") * attempt
	shell:
		"generate_metrics {input[0]} > {output[0]} 2> {log}"

rule minos_mikado_serialise:
	input:
		config = config["mikado-config-file"],
		ext_scores = rules.minos_metrics_generate_metrics_matrix.output[0],
		transcripts = rules.minos_mikado_prepare.output[0]
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
		HPC_CONFIG.get_cores("minos_mikado_serialise")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_mikado_serialise") * attempt
	shell:
		"{params.program_call} {params.program_params} --transcripts {input.transcripts} --external-scores {input.ext_scores} --json-conf {input.config} --procs {threads} -od {params.outdir} &> {log}" + \
		" && touch {output[0]}"

rule minos_mikado_pick:
	input:
		config = config["mikado-config-file"],
		gtf = rules.minos_mikado_prepare.output[1],
		serialise_done = rules.minos_mikado_serialise.output[0],
		db = rules.minos_mikado_serialise.output[1]
	output:
		loci = os.path.join(config["outdir"], "mikado.loci.gff3"),
		subloci = os.path.join(config["outdir"], "mikado.subloci.gff3"),
		monoloci = os.path.join(config["outdir"], "mikado.monoloci.gff3"),
		loci_metrics = os.path.join(config["outdir"], "mikado.loci.metrics.tsv"),
		subloci_metrics = os.path.join(config["outdir"], "mikado.subloci.metrics.tsv"),
		monoloci_metrics = os.path.join(config["outdir"], "mikado.monoloci.metrics.tsv")
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="pick"),
		program_params = config["params"]["mikado"]["pick"],
		outdir = config["outdir"]
	threads:
		HPC_CONFIG.get_cores("minos_mikado_pick")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_mikado_pick") * attempt
	shell:
		"{params.program_call} {params.program_params}" + \
		" -od {params.outdir} --procs {threads} --json-conf {input.config}" + \
		" --subloci-out $(basename {output.subloci})" + \
		" --monoloci-out $(basename {output.monoloci})" + \
		" -db {input.db} {input.gtf}"

rule minos_parse_mikado_pick:
	input:
		loci = rules.minos_mikado_pick.output[0]
	output:
		gff = os.path.join(config["outdir"], POST_PICK_PREFIX + ".gff")
	shell:
		"parse_mikado_gff {input.loci} > {output.gff}"

rule minos_gffread_extract_sequences_post_pick:
	input:
		gff = rules.minos_parse_mikado_pick.output[0],
		refseq = config["reference-sequence"]
	output:
		cdna = os.path.join(config["outdir"], POST_PICK_PREFIX + ".cdna.fasta"),
		tbl = os.path.join(config["outdir"], POST_PICK_PREFIX + ".table.txt"),
		cds = os.path.join(config["outdir"], POST_PICK_PREFIX + ".cds.fasta"),
		pep = os.path.join(config["outdir"], POST_PICK_PREFIX + ".proteins.fasta"),
	params:
		program_call = config["program_calls"]["gffread"],
		table_format = "--table @chr,@start,@end,@strand,@numexons,@covlen,@cdslen,ID,Note,confidence,representative,biotype,InFrameStop,partialness"
	shell:
		"{params.program_call} {input.gff} -g {input.refseq} -P {params.table_format} -W -w {output.cdna} -x {output.cds} -y {output.pep} -o {output.tbl}"


rule minos_calculate_cds_lengths_post_pick:
	input:
		cds = rules.minos_gffread_extract_sequences_post_pick.output.cds,
		cdna = rules.minos_gffread_extract_sequences_post_pick.output.cdna
	output:
		rules.minos_gffread_extract_sequences_post_pick.output.cds + ".lengths"
	params:
		min_cds_length = config["misc"]["min_cds_length"]
	run:
		from minos.scripts.calculate_cdslen import calculate_cdslen
		calculate_cdslen(input.cds, input.cdna, output[0], params.min_cds_length)


rule minos_gff_genometools_check_post_pick:
	input:
		rules.minos_parse_mikado_pick.output[0]
	output:
		gff = os.path.join(config["outdir"], POST_PICK_PREFIX + ".gt_checked.gff")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".post_pick_genometools_check.log")
	params:
		program_call = config["program_calls"]["genometools"],
		program_params = config["params"]["genometools"]["check"]
	shell:
		"{params.program_call} {params.program_params} {input[0]} > {output.gff} 2> {log}"

rule minos_gff_validate_post_gt:
	input:
		rules.minos_gff_genometools_check_post_pick.output[0]
	output:
		os.path.join(config["outdir"], POST_PICK_PREFIX + ".gt_checked.validation_report.txt")
	threads:
		HPC_CONFIG.get_cores("minos_gff_validate_post_gt")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_gff_validate_post_gt") * attempt
	shell:
		"validate_gff3 {input} > {output}"

rule minos_kallisto_index_post_pick:
	input:
		rules.minos_gffread_extract_sequences_post_pick.output.cdna
	output:
		os.path.join(config["outdir"], "kallisto", os.path.basename(rules.minos_gffread_extract_sequences_post_pick.output.cdna) + ".idx")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".kallisto_index_post_pick.log")
	params:
		program_call = config["program_calls"]["kallisto"].format(program="index"),
		program_params = config["params"].get("kallisto", {}).get("index", "")
	threads:
		HPC_CONFIG.get_cores("minos_kallisto_index_post_pick")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_kallisto_index_post_pick") * attempt
	shell:
		"{params.program_call} {params.program_params} -i {output[0]} {input[0]} &> {log}"

rule minos_kallisto_quant_post_pick:
	input:
		index = rules.minos_kallisto_index_post_pick.output[0],
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
		HPC_CONFIG.get_cores("minos_kallisto_quant_post_pick")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_kallisto_quant_post_pick") * attempt
	shell:
		"{params.program_call} {params.program_params} {params.stranded} -i {input.index} -o {params.outdir} --threads {threads} {input.reads} &> {log}"

rule minos_collapse_metrics:
	input:
		gff = rules.minos_parse_mikado_pick.output[0],
		ext_scores = rules.minos_metrics_generate_metrics_matrix.output[0],
		metrics_info = rules.minos_metrics_generate_metrics_info.output[0],
		expression = expand(rules.minos_kallisto_quant_post_pick.output, run=config["data"]["expression-runs"].keys()),
		cds_lengths = rules.minos_calculate_cds_lengths_post_pick.output[0]
	output:
		os.path.join(config["outdir"], POST_PICK_PREFIX + ".collapsed_metrics.tsv"),
		os.path.join(config["outdir"], "COLLAPSE_METRICS_DONE"),
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".collapse_metrics.log")
	threads:
		HPC_CONFIG.get_cores("minos_collapse_metrics")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_collapse_metrics") * attempt
	run:
		from minos.scripts.collapse_metrics import MetricCollapser
		mc = MetricCollapser(input.gff, input.metrics_info, input.ext_scores, input.cds_lengths, input.expression)
		with open(output[0], "w") as out:
			mc.write_scores(config["collapse_metrics_thresholds"], stream=out)
		open(output[1], "w").close()

rule minos_summarise_collapsed_metrics:
	input:
		rules.minos_collapse_metrics.output[0]
	output:
		os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".release.collapsed_metrics.summary.tsv")
	run:
		from minos.scripts.summarise_collapsed_metrics import CollapsedMetricsSummariser
		with open(output[0], "w") as out:
			CollapsedMetricsSummariser(input[0]).write_summary(stream=out)


rule minos_create_release_gffs:
	input:
		gff = rules.minos_gff_genometools_check_post_pick.output[0],
		metrics_info = rules.minos_collapse_metrics.output[0],
		sentinel = rules.minos_collapse_metrics.output[1]
	output:
		os.path.join(config["outdir"], POST_PICK_PREFIX + ".release.unsorted.gff3"),
		os.path.join(config["outdir"], POST_PICK_PREFIX + ".release_browser.unsorted.gff3"),
		rules.minos_gff_genometools_check_post_pick.output[0] + ".old_new_id_relation.txt"
	params:
		annotation_version = config.get("annotation_version", "EIv1"),
		genus_identifier = config.get("genus_identifier", "XYZ")
	threads:
		HPC_CONFIG.get_cores("minos_create_release_gffs")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_create_release_gffs") * attempt
	shell:
		"create_release_gff3 {input.gff} {input.metrics_info} --annotation-version {params.annotation_version} --genus-identifier {params.genus_identifier} 2> {LOG_DIR}/create_release_gff.log"

rule minos_create_release_metrics:
	input:
		rules.minos_create_release_gffs.output[2],
		rules.minos_collapse_metrics.output[0]
	output:
		os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".release.metrics.tsv")
	run:
		import csv
		genes, transcripts = dict(), dict()
		for new_gene, new_transcript, old_gene, old_transcript in csv.reader(open(input[0]), delimiter="\t"):
			if not new_gene.startswith("#"):
				genes[old_gene] = new_gene
				transcripts[old_transcript] = new_transcript
		with open(output[0], "w") as metrics_out:
			for row in csv.reader(open(input[1]), delimiter="\t"):
				if not row[0].startswith("#"):
					row[1] = genes.get(row[1])
					row[0] = transcripts.get(row[0])
					if row[0] is None:
						continue
				print(*row, sep="\t", flush=True, file=metrics_out)

rule minos_sort_release_gffs:
	input:
		rules.minos_create_release_gffs.output
	output:
		os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".release.gff3"),
		os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".release.browser.gff3")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".sort_release_gffs.log")
	params:
		program_call = config["program_calls"]["genometools"],
		program_params = config["params"]["genometools"]["sort"]
	threads:
		HPC_CONFIG.get_cores("minos_sort_release_gffs")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_sort_release_gffs") * attempt
	shell:
		"{params.program_call} {params.program_params} {input[0]} > {output[0]} 2> {log}" + \
		" && {params.program_call} {params.program_params} {input[1]} > {output[1]} 2>> {log}"

rule minos_final_sanity_check:
	input:
		rules.minos_sort_release_gffs.output[0]
	output:
		temp(os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".sanity_checked.release.gff3"))
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".final_sanity_check.log")
	threads:
		HPC_CONFIG.get_cores("minos_final_sanity_check")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_final_sanity_check") * attempt
	shell:
		"sanity_check {input[0]} > {output[0]} 2> {log}"

rule minos_mikado_pick_extract_coords:
	input:
		rules.minos_sort_release_gffs.output[0]
	output:
		rules.minos_sort_release_gffs.output[0] + ".coords"
	run:
		from minos.scripts.extract_coords import extract_coords
		extract_coords(input[0], output[0], filetype="gff")

rule minos_generate_mikado_stats:
	input:
		rules.minos_final_sanity_check.output[0]
	output:
		rules.minos_sort_release_gffs.output[0] + ".mikado_stats.txt",
		rules.minos_sort_release_gffs.output[0] + ".mikado_stats.tsv"
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="util stats"),
	threads:
		HPC_CONFIG.get_cores("minos_generate_mikado_stats")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_generate_mikado_stats") * attempt
	shell:
		"{params.program_call} {input} --tab-stats {output[1]} > {output[0]}" + \
		" && parse_mikado_stats {output[0]} > {output[0]}.summary"

rule minos_extract_final_sequences:
	input:
		gff = rules.minos_final_sanity_check.output[0],
		refseq = config["reference-sequence"]
	output:
		cdna = rules.minos_sort_release_gffs.output[0] + ".cdna.fasta",
		tbl = rules.minos_sort_release_gffs.output[0] + ".gffread.table.txt",
		cds = rules.minos_sort_release_gffs.output[0] + ".cds.fasta",
		pep = rules.minos_sort_release_gffs.output[0] + ".pep.raw.fasta"
	params:
		program_call = config["program_calls"]["gffread"],
		table_format = "--table @chr,@start,@end,@strand,@numexons,@covlen,@cdslen,ID,Note,confidence,representative,biotype,InFrameStop,partialness"
	shell:
		"{params.program_call} {input.gff} -g {input.refseq} -P {params.table_format} -W -w {output.cdna} -x {output.cds} -y {output.pep} -o {output.tbl}"

rule minos_cleanup_final_proteins:
	input:
		rules.minos_extract_final_sequences.output.pep
	output:
		rules.minos_extract_final_sequences.output.pep.replace(".raw.fasta", ".fasta")
	log:
		os.path.join(LOG_DIR, "cleanup_proteins.log")
	params:
		prefix = rules.minos_extract_final_sequences.output.pep.replace(".raw.fasta", ""),
		program_call = config["program_calls"]["prinseq"],
		program_params = config["params"]["prinseq"]
	threads:
		HPC_CONFIG.get_cores("minos_cleanup_final_proteins")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_cleanup_final_proteins") * attempt
	shell:
		"{params.program_call} -aa -fasta {input} {params.program_params} -out_good {params.prefix} -out_bad {params.prefix}.bad"

rule minos_generate_final_table:
	input:
		stats_table = rules.minos_generate_mikado_stats.output[1],
		seq_table = rules.minos_extract_final_sequences.output.tbl,
		bt_conf_table = rules.minos_final_sanity_check.output[0]
	output:
		final_table = rules.minos_sort_release_gffs.output[0] + ".final_table.tsv",
		summary = rules.minos_sort_release_gffs.output[0] + ".biotype_conf.summary"
	threads:
		HPC_CONFIG.get_cores("minos_generate_final_table")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("minos_generate_final_table") * attempt
	run:
		from minos.scripts.generate_final_table import generate_final_table
		generate_final_table(input.seq_table, input.bt_conf_table, input.stats_table, output.final_table, output.summary)

rule minos_collate_metric_oddities:
	input:
		loci = rules.minos_mikado_pick.output.loci_metrics,
		subloci = rules.minos_mikado_pick.output.subloci_metrics,
		monoloci = rules.minos_mikado_pick.output.monoloci_metrics,
		old_new_rel = rules.minos_create_release_gffs.output[2],
		final_table = rules.minos_generate_final_table.output.final_table,
	output:
		os.path.join(config["outdir"], "results", RELEASE_PREFIX + ".release.metrics_oddities.tsv"),
	run:
		import csv
		from minos.scripts.metric_oddities import MetricOddityParser
		release_transcripts = {row[0]: row[14:16] for row in csv.reader(open(input.final_table), delimiter="\t") if not row[0].startswith("#")}
		transcript_data = {
			row[3]: release_transcripts.get(row[1])
			for row in csv.reader(open(input.old_new_rel), delimiter="\t") if not row[0].startswith("#") and row[1] in release_transcripts
		}
		with open(output[0], "w") as loci_oddities_out:
			MetricOddityParser(input[0], input[1], input[2], config["report_metric_oddities"], transcript_data=transcript_data).write_table(stream=loci_oddities_out)

rule config_copy_results:
	input:
		rules.minos_collate_metric_oddities.output[0],
		CONFIG_COPY_SOURCES
	output:
		CONFIG_COPY_TARGETS
	run:
		import pathlib
		import os
		import shutil
		for src, tgt in CONFIG_COPY:
			tgt_dir = os.path.dirname(tgt)
			pathlib.Path(tgt_dir).mkdir(exist_ok=True, parents=True)
			shutil.copyfile(src, tgt)

rule split_proteins_prepare:
	input:
		rules.minos_gffread_extract_sequences.output[0]
	output:
		expand(os.path.join(BUSCO_PATH, "runs", "proteins_prepare", "input", "{run}.proteins.fasta"), run=config["data"]["transcript_models"])
	log:
		os.path.join(BUSCO_PATH, "logs", "split_proteins_prepare.log")
	run:
		from minos.scripts.busco_splitter import split_fasta
		fasta_files = {tm: open(os.path.join(BUSCO_PATH, "runs", "proteins_prepare", "input", tm + ".proteins.fasta"), "w") for tm in config["data"]["transcript_models"]}
		split_fasta(input[0], fasta_files)

rule busco_proteins_prepare:
	input:
		os.path.join(BUSCO_PATH, "runs", "proteins_prepare", "input", "{run}.proteins.fasta")
	output:
		os.path.join(BUSCO_PATH, "runs", "proteins_prepare", "{run}", "run_{}".format(BUSCO_LINEAGE), "short_summary.txt"),
		os.path.join(BUSCO_PATH, "runs", "proteins_prepare", "{run}", "run_{}".format(BUSCO_LINEAGE), "full_table.tsv"),
		os.path.join(BUSCO_PATH, "runs", "proteins_prepare", "{run}", "run_{}".format(BUSCO_LINEAGE), "missing_busco_list.tsv")
	log:
		os.path.join(BUSCO_PATH, "logs", "{run}.proteins_prepare.log")
	params:
		input = lambda wildcards: os.path.join(BUSCO_PATH, "runs", "proteins_prepare", "input", wildcards.run + ".proteins.fasta"),
		program_call = config["program_calls"]["busco"],
		program_params = config["params"]["busco"]["proteins_prepare"],
		lineage_path = config["busco_analyses"]["lineage"],
		run = lambda wildcards: wildcards.run,
		busco_mode = "proteins",
		busco_stage = "proteins_prepare",
		copy_busco = config["program_calls"]["copy_busco"]
	threads:
		HPC_CONFIG.get_cores("busco_proteins_prepare")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("busco_proteins_prepare") * attempt
	shell:
		BUSCO_CMD

rule busco_proteins_final:
	input:
		rules.minos_extract_final_sequences.output.pep
	output:
		os.path.join(BUSCO_PATH, "runs", "proteins_final", "proteins_final", "run_{}".format(BUSCO_LINEAGE), "short_summary.txt"),
		os.path.join(BUSCO_PATH, "runs", "proteins_final", "proteins_final", "run_{}".format(BUSCO_LINEAGE), "full_table.tsv"),
		os.path.join(BUSCO_PATH, "runs", "proteins_final", "proteins_final", "run_{}".format(BUSCO_LINEAGE), "missing_busco_list.tsv")
	log:
		os.path.join(BUSCO_PATH, "logs", "proteins_final.log")
	params:
		input = os.path.abspath(rules.minos_extract_final_sequences.output.pep),
		program_call = config["program_calls"]["busco"],
		program_params = config["params"]["busco"]["proteins_final"],
		lineage_path = config["busco_analyses"]["lineage"],
		run = "proteins_final",
		busco_mode = "proteins",
		busco_stage = "proteins_final",
		copy_busco = config["program_calls"]["copy_busco"]
	threads:
		HPC_CONFIG.get_cores("busco_proteins_final")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("busco_proteins_final") * attempt
	shell:
		BUSCO_CMD


rule split_transcripts_prepare:
	input:
		rules.minos_gffread_extract_sequences.output[1]
	output:
		expand(os.path.join(BUSCO_PATH, "runs", "transcripts_prepare", "input", "{run}.cdna.fasta"), run=config["data"]["transcript_models"])
	log:
		os.path.join(BUSCO_PATH, "logs", "split_transcripts_prepare.log")
	run:
		from minos.scripts.busco_splitter import split_fasta
		fasta_files = {tm: open(os.path.join(BUSCO_PATH, "runs", "transcripts_prepare", "input", tm + ".cdna.fasta"), "w") for tm in config["data"]["transcript_models"]}
		split_fasta(input[0], fasta_files)

rule busco_transcripts_prepare:
	input:
		os.path.join(BUSCO_PATH, "runs", "transcripts_prepare", "input", "{run}.cdna.fasta")
	output:
		os.path.join(BUSCO_PATH, "runs", "transcripts_prepare", "{run}", "run_{}".format(BUSCO_LINEAGE), "short_summary.txt"),
		os.path.join(BUSCO_PATH, "runs", "transcripts_prepare", "{run}", "run_{}".format(BUSCO_LINEAGE), "full_table.tsv"),
		os.path.join(BUSCO_PATH, "runs", "transcripts_prepare", "{run}", "run_{}".format(BUSCO_LINEAGE), "missing_busco_list.tsv")
	log:
		os.path.join(BUSCO_PATH, "logs", "{run}.transcripts_prepare.log")
	params:
		input = lambda wildcards: os.path.join(BUSCO_PATH, "runs", "transcripts_prepare", "input", wildcards.run + ".cdna.fasta"),
		program_call = config["program_calls"]["busco"],
		program_params = config["params"]["busco"]["transcripts_prepare"],
		lineage_path = config["busco_analyses"]["lineage"],
		run = lambda wildcards: wildcards.run,
		busco_mode = "transcriptome",
		busco_stage = "transcripts_prepare",
		copy_busco = config["program_calls"]["copy_busco"]
	threads:
		HPC_CONFIG.get_cores("busco_transcripts_prepare")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("busco_transcripts_prepare") * attempt
	shell:
		BUSCO_CMD

rule busco_transcripts_final:
	input:
		rules.minos_extract_final_sequences.output.cdna
	output:
		os.path.join(BUSCO_PATH, "runs", "transcripts_final", "transcripts_final", "run_{}".format(BUSCO_LINEAGE), "short_summary.txt"),
		os.path.join(BUSCO_PATH, "runs", "transcripts_final", "transcripts_final", "run_{}".format(BUSCO_LINEAGE), "full_table.tsv"),
		os.path.join(BUSCO_PATH, "runs", "transcripts_final", "transcripts_final", "run_{}".format(BUSCO_LINEAGE), "missing_busco_list.tsv"),
	log:
		os.path.join(BUSCO_PATH, "logs", "transcripts_final.log")
	params:
		input = os.path.abspath(rules.minos_extract_final_sequences.output.cdna),
		program_call = config["program_calls"]["busco"],
		program_params = config["params"]["busco"]["transcripts_final"],
		lineage_path = config["busco_analyses"]["lineage"],
		run = "transcripts_final",
		busco_mode = "transcriptome",
		busco_stage = "transcripts_final",
		copy_busco = config["program_calls"]["copy_busco"]
	threads:
		HPC_CONFIG.get_cores("busco_transcripts_final")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("busco_transcripts_final") * attempt
	shell:
		BUSCO_CMD

if config["busco_analyses"]["genome"] or config["busco_analyses"]["precomputed_genome"]:
	rule busco_genome:
		input:
			config["busco_analyses"]["precomputed_genome"]["summary"] if config["busco_analyses"]["precomputed_genome"] else config["reference-sequence"],
			config["busco_analyses"]["precomputed_genome"]["full_table"] if config["busco_analyses"]["precomputed_genome"] else config["reference-sequence"],
			config["busco_analyses"]["precomputed_genome"]["missing_busco_list"] if config["busco_analyses"]["precomputed_genome"] else config["reference-sequence"]
		output:
			BUSCO_GENOME_OUTPUT_SUMMARY,
			BUSCO_GENOME_OUTPUT_FULL_TABLE,
			BUSCO_GENOME_MISSING
		log:
			os.path.join(BUSCO_PATH, "logs", "genome.log")
		params:
			input = os.path.abspath(config["reference-sequence"]),
			program_call = config["program_calls"]["busco"],
			program_params = config["params"]["busco"]["genome"],
			lineage_path = config["busco_analyses"]["lineage"],
			run = "genome",
			busco_mode = "genome",
			busco_stage = "genome",
			copy_busco = config["program_calls"]["copy_busco"]
		threads:
			HPC_CONFIG.get_cores("busco_genome")
		resources:
			mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("busco_genome") * attempt
		shell:
			BUSCO_PRECOMPUTED if config["busco_analyses"]["precomputed_genome"] else BUSCO_CMD

rule busco_copy_results:
	input:
		BUSCO_COPY_SOURCES
	output:
		BUSCO_COPY_TARGETS
	run:
		import pathlib
		import os
		import shutil
		for src, tgt in BUSCO_COPY:
			tgt_dir = os.path.dirname(tgt)
			pathlib.Path(tgt_dir).mkdir(exist_ok=True, parents=True)
			shutil.copyfile(src, tgt)

rule busco_summary:
	input:
		rules.minos_mikado_prepare_extract_coords.output[0],
		rules.minos_mikado_pick_extract_coords.output[0],
		BUSCO_ANALYSES + BUSCO_PROTEIN_PREPARE_RUNS
	output:
		BUSCO_TABLE
	run:
		from minos.scripts.generate_busco_tables import BuscoTableGenerator

		btg = BuscoTableGenerator(
			os.path.join(config["outdir"], "tx2gene"),
			input[0],
			input[1],
			os.path.join(BUSCO_PATH, "runs")
		)
		btg.write_review_table(output[0])
		btg.write_raw_data(output[0])
		btg.write_busco_table(output[0], config["misc"]["busco_max_copy_number"])
