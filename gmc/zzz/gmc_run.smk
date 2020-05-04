import os
import sys
import pathlib

from gmc.gmc_configure import ExternalMetrics
from eicore.hpc_config import HpcConfig
HPC_CONFIG = HpcConfig(config["hpc_config"])

EXTERNAL_METRICS_DIR = os.path.join(config["outdir"], "generate_metrics")
pathlib.Path(EXTERNAL_METRICS_DIR).mkdir(exist_ok=True, parents=True)

LOG_DIR = os.path.join(config["outdir"], "logs")
TEMP_DIR = os.path.join(config["outdir"], "tmp")
AUGUSTUS_CONFIG_DATA = config["paths"]["augustus_config_data"]


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
	return config["transcript_models"][wc.run]

# output management
POST_PICK_PREFIX = "mikado.annotation"
RELEASE_PREFIX = config.get("genus_identifier", "XYZ") + "_" + config.get("annotation_version", "EIv1")
RESULTS_DIR = os.path.join(config["outdir"], "results")

OUTPUTS = [
	os.path.join(config["outdir"], "mikado_prepared.fasta"),
	os.path.join(config["outdir"], "mikado_prepared.exon.gff"),
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

for run in config.get("data", dict()).get("protein-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "proteins", run, "mikado_" + run + ".refmap"))
for run in config.get("data", dict()).get("transcript-runs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "mikado_compare", "transcripts", run, "mikado_" + run + ".refmap"))
for run in config.get("data", dict()).get("protein-seqs", dict()):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], run, run + ".{}.tsv.tophit".format(config["blast-mode"])))
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
	&& {params.copy} {AUGUSTUS_CONFIG_DATA} $cfgdir/
	&& export AUGUSTUS_CONFIG_PATH=$cfgdir/config
	&& rm -rvf {BUSCO_PATH}/runs/{params.busco_stage}/{params.run}
	&& mkdir -p {BUSCO_PATH}/runs/{params.busco_stage}/{params.run}
	&& cd {BUSCO_PATH}/runs/{params.busco_stage}/{params.run}
	&& {params.program_call} {params.program_params} -i {params.input} -c {threads} -m {params.busco_mode} --force -l {params.lineage_path} -o {params.run} &> {log}
	&& mv {params.run}/* . && rm -rf {params.run}
	&& touch {output[2]}
	&& rm -rf $cfgdir
""".strip().replace("\n\t", " ")

TX2GENE_MAPS = [
	os.path.join(config["outdir"], "tx2gene", tm + ".tx2gene")
	for tm in config["transcript_models"]
]

BUSCO_PATH = os.path.abspath(os.path.join(config["outdir"], "busco"))
BUSCO_LINEAGE = os.path.basename(config["busco_analyses"]["lineage"]) if config["busco_analyses"]["lineage"] is not None else None

BUSCO_ANALYSES = list()
BUSCO_PROTEIN_PREPARE_RUNS = list()
if config["busco_analyses"]["proteins"]:
	BUSCO_PROTEIN_PREPARE_RUNS.extend(
		os.path.join(BUSCO_PATH, "runs", "proteins_prepare", "{tm}", "run_{lineage}", busco_file).format(tm=tm, lineage=BUSCO_LINEAGE)
		for tm in config["transcript_models"] for busco_file in ("short_summary.txt", "full_table.tsv", "missing_busco_list.tsv")
	)
	BUSCO_ANALYSES.extend(
		os.path.join(BUSCO_PATH, "runs", "proteins_final", "proteins_final", "run_{lineage}", busco_file).format(lineage=BUSCO_LINEAGE)
		for busco_file in ("short_summary.txt", "full_table.tsv", "missing_busco_list.tsv")
	)

if config["busco_analyses"]["transcriptome"]:
	BUSCO_ANALYSES.extend(
		os.path.join(BUSCO_PATH, "runs", "transcripts_prepare", "{tm}", "run_{lineage}", busco_file).format(tm=tm, lineage=BUSCO_LINEAGE)
		for tm in config["transcript_models"] for busco_file in ("short_summary.txt", "full_table.tsv", "missing_busco_list.tsv")
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

if BUSCO_ANALYSES or BUSCO_PROTEIN_PREPARE_RUNS:
	BUSCO_ANALYSES.extend(TX2GENE_MAPS)
	BUSCO_TABLE = os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".busco_final_table.tsv")
	OUTPUTS.extend(BUSCO_PROTEIN_PREPARE_RUNS)

localrules:
	all,
	gmc_extract_exons,
	gmc_metrics_repeats_convert,
	gmc_metrics_repeats_extract_exons,
	gmc_mikado_prepare_extract_coords,
	gmc_mikado_pick_extract_coords,
	gmc_metrics_parse_repeat_coverage,
	gmc_metrics_blastp_combine,
	gmc_metrics_generate_metrics_info,
	gmc_parse_mikado_pick,
	gmc_gffread_extract_sequences_post_pick,
	gmc_gffread_extract_sequences,
	gmc_gff_genometools_check_post_pick,
	gmc_calculate_cds_lengths_post_pick,
	gmc_extract_final_sequences,
	split_proteins_prepare,
	split_transcripts_prepare,
	busco_copy_results,
	busco_concat_protein_metrics,
	busco_summary


rule all:
	input:
		OUTPUTS,
		os.path.join(EXTERNAL_METRICS_DIR, "metrics_info.txt"),
		os.path.join(config["outdir"], "MIKADO_SERIALISE_DONE"),
		os.path.join(config["outdir"], "mikado.subloci.gff3"),
		os.path.join(config["outdir"], "mikado.loci.gff3"),
		[
			os.path.join(config["outdir"], POST_PICK_PREFIX + suffix)
			for suffix in {".gff", ".proteins.fasta", ".table.txt", ".cds.fasta", ".cds.fasta.lengths", ".cdna.fasta"}
		],
		os.path.join(config["outdir"], "kallisto", POST_PICK_PREFIX + ".cdna.fasta.idx"),
		POST_PICK_EXPRESSION,
		[
			os.path.join(config["outdir"], POST_PICK_PREFIX + suffix)
			for suffix in {".gt_checked.gff", ".collapsed_metrics.tsv", ".gt_checked.validation_report.txt", ".release.unsorted.gff3", ".release_browser.unsorted.gff3"}
		],
		[
			os.path.join(RESULTS_DIR, RELEASE_PREFIX + suffix)
			for suffix in {".release.gff3", ".release_browser.gff3",
							".sanity_checked.release.gff3", ".sanity_checked.release.gff3.mikado_stats.txt", ".sanity_checked.release.gff3.mikado_stats.tsv"}
		],
		[
			os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".sanity_checked.release.gff3") + ".{}.fasta".format(dtype) for dtype in {"cdna", "cds", "pep.raw", "pep"}
		],
		[
			os.path.join(RESULTS_DIR, RELEASE_PREFIX + suffix)
			for suffix in {".sanity_checked.release.gff3.final_table.tsv", ".sanity_checked.release.gff3.biotype_conf.summary"}
		],
		BUSCO_ANALYSES + ([BUSCO_TABLE] if BUSCO_TABLE else []),
		BUSCO_COPY_TARGETS

rule gmc_mikado_prepare:
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
		HPC_CONFIG.get_cores("gmc_mikado_prepare")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_mikado_prepare") * attempt
	shell:
		"{params.program_call} {params.program_params} --json-conf {input[0]} --procs {threads} -od {params.outdir} &> {log}"


rule gmc_mikado_prepare_extract_coords:
	input:
		rules.gmc_mikado_prepare.output[1]
	output:
		rules.gmc_mikado_prepare.output[1] + ".coords"
	run:
		import csv
		import re
		with open(output[0], "w") as coords_out:
			for row in csv.reader(open(input[0]), delimiter="\t"):
				if row and not row[0].startswith("#") and "rna" in row[2].lower():
					tid = re.search('transcript_id "?([^";]+)"?;', row[8]).group(1)
					print(tid, "{seq}:{start}..{end}".format(seq=row[0], start=row[3], end=row[4]), sep="\t", flush=True, file=coords_out)

rule gmc_generate_tx2gene_maps:
	input:
		get_transcript_models
	output:
		os.path.join(config["outdir"], "tx2gene", os.path.basename("{run}") + ".tx2gene")
	threads:
		HPC_CONFIG.get_cores("gmc_generate_tx2gene_maps")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_generate_tx2gene_maps")
	run:
		import csv
		import re
		with open(output[0], "w") as tx2gene_out:
			for row in csv.reader(open(input[0]), delimiter="\t"):
				if not row or (row and row[0].startswith("#")):
					continue
				if row[2] in {"mRNA", "ncRNA"}:
					attr = dict((item.group(1).strip(), item.group(2).strip()) for item in re.finditer("([^;]+)\s*=\s*([^;]+);?", row[8]))
					print("{}_{}".format(wildcards.run, attr["ID"]), attr["Parent"], file=tx2gene_out, flush=True, sep="\t")

rule gmc_extract_exons:
	input:
		rules.gmc_mikado_prepare.output[1]
	output:
		rules.gmc_mikado_prepare.output[1].replace(".gtf", ".exon.gff")
	shell:
		"""
		awk '$3 == "exon"' {input[0]} | awk -v FS="\\t" -v OFS="\\t" '{{exon += 1; split($9,a,"\\""); $9="ID="a[4]".exon"exon";Parent="a[4]; print $0}}' > {output[0]}
		""".strip().replace("\n\t", " ")

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
	threads:
		HPC_CONFIG.get_cores("gmc_mikado_compare_index_reference")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_mikado_compare_index_reference") * attempt
	shell:
		"{params.program_call} {params.program_params} -r {input} &> {log}"

rule gmc_gffread_extract_sequences:
	input:
		gtf = rules.gmc_mikado_prepare.output[1],
		refseq = config["reference-sequence"]
	output:
		rules.gmc_mikado_prepare.output[1] + (".prot.fasta" if config["blast-mode"] == "blastp" else ".cds.fasta"),
		rules.gmc_mikado_prepare.output[1] + ".cdna.fasta"
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".gffread_extract.log")
	params:
		program_call = config["program_calls"]["gffread"],
		program_params = config["params"]["gffread"][config["blast-mode"]],
		output_params = "-W -x" if config["blast-mode"] == "blastx" else ("-y" if config["blast-mode"] == "blastp" else "")
	shell:
		"{params.program_call} {input.gtf} -g {input.refseq} {params.program_params} -W -w {output[1]} {params.output_params} {output[0]} &> {log}"

rule gmc_metrics_repeats_convert:
	input:
		get_repeat_data
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "repeats", "{run}.converted.gff")
	run:
		import csv
		import re
		source = tag = wildcards.run
		name_tag, rep_counter = "RM", 0
		with open(output[0], "w") as out_gff:
			for row in csv.reader(open(input[0]), delimiter="\t"):
				if not row or (row and row[0].startswith("#")):
					continue
				if row[1] == "RepeatMasker":
					rep_counter += 1
					try:
						target = re.search('Target\s+"([^"]+)', row[8]).group(1)
					except:
						try:
				 			target = re.search('Target\s*=\s*([^;]+)', row[8]).group(1)
						except:
							raise ValueError("Cannot parse Target from " + row[8])
					target = re.sub("\s+", "_", target)
					note = ";Note={}".format(target)
					try:
						name = re.sub("\s+", "_", re.search('Name\s*=\s*([^;]+)', row[8]).group(1))
					except:
						name, note = target, ""
					row[1] = source
					row[5] = "{:0.0f}".format(float(row[5]))
					row[2] = "match"
					attrib = "ID={tag}:{name_tag}{counter};Name={name}{note}".format(tag=tag, name_tag=name_tag, counter=rep_counter, name=name, note=note)
					print(*row[:8], attrib, sep="\t", file=out_gff, flush=True)	
					row[2] = "match_part"
					attrib = "ID={tag}:{name_tag}{counter}-exon1;Parent={tag}:{name_tag}{counter}".format(tag=tag, name_tag=name_tag, counter=rep_counter)
					print(*row[:8], attrib, sep="\t", file=out_gff, flush=True)	
					print("###", file=out_gff, flush=True)
			

rule gmc_metrics_repeats_extract_exons:
	input:
		rules.gmc_metrics_repeats_convert.output[0]
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "repeats", "{run}.no_strand.exon.gff")
	shell:
		"""
		awk '$3 == "match_part"' {input[0]} | sed -e 's/\\tmatch_part\\t/\\texon\\t/' -e 's/\\t[+-]\\t/\\t.\\t/' > {output[0]}
		""".strip().replace("\n\t", " ")

rule gmc_metrics_bedtools_repeat_coverage:
	input:
		rules.gmc_metrics_repeats_extract_exons.output[0],
		rules.gmc_extract_exons.output[0]
	output:
		rules.gmc_metrics_repeats_extract_exons.output[0] + ".cbed",
	params:
		program_call = config["program_calls"]["bedtools"]["coverageBed"]
	threads:
		HPC_CONFIG.get_cores("gmc_metrics_bedtools_repeat_coverage")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_metrics_bedtools_repeat_coverage") * attempt
	shell:
		"""
		{params.program_call} -a {input[1]} -b {input[0]} > {output[0]}.tmp
		&& cut -f 9 {output[0]}.tmp | cut -d ';' -f 1 | paste - {output[0]}.tmp | sort -k1,1V | cut -f 2- > {output[0]}
		&& rm {output[0]}.tmp
		""".strip().replace("\n\t", " ")


rule gmc_metrics_parse_repeat_coverage:
	input:
		rules.gmc_metrics_bedtools_repeat_coverage.output[0],
	output:
		rules.gmc_metrics_bedtools_repeat_coverage.output[0] + ".parsed.txt",
	run:
		from gmc.scripts.parse_cbed_stats import parse_cbed
		with open(output[0], "w") as outstream, open(input[0]) as instream:
			parse_cbed(instream, outstream=outstream)


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
		HPC_CONFIG.get_cores("gmc_metrics_cpc2")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_metrics_cpc2") * attempt	
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
		threads:
			HPC_CONFIG.get_cores("gmc_metrics_kallisto_index")
		resources:
			mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_metrics_kallisto_index") * attempt
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
			HPC_CONFIG.get_cores("gmc_metrics_kallisto_quant")
		resources:
			mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_metrics_kallisto_quant") * attempt
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
	threads:
		HPC_CONFIG.get_cores("gmc_metrics_mikado_compare_vs_transcripts")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_metrics_mikado_compare_vs_transcripts") * attempt
	shell:
		"mkdir -p {params.outdir}" + \
		" && {params.program_call} {params.program_params} -r {input.mika} -p {input.transcripts} -o {params.outdir}/mikado_{params.transcripts} &> {log}" + \
		" && touch {output[0]}"
		
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
	threads:
		HPC_CONFIG.get_cores("gmc_metrics_mikado_compare_vs_proteins")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_metrics_mikado_compare_vs_proteins") * attempt
	shell:
		"mkdir -p {params.outdir}" + \
		" && {params.program_call} {params.program_params} -r {input.mika} -p {input.proteins} -o {params.outdir}/mikado_{params.proteins} &> {log}" + \
		" && touch {output[0]}"

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
	threads:
		HPC_CONFIG.get_cores("gmc_metrics_blastp_mkdb")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_metrics_blastp_mkdb") * attempt
	shell:
		"{params.program_call} {params.program_params} -in {input[0]} -out {params.outdir}/{params.db_prefix} -logfile {log}" + \
		" && touch {output[0]}"

checkpoint gmc_chunk_proteins:
	input:
		rules.gmc_gffread_extract_sequences.output[0]
	output:
		chunk_dir = directory(os.path.join(TEMP_DIR, "chunked_proteins"))
	log:
		os.path.join(LOG_DIR, os.path.basename(rules.gmc_gffread_extract_sequences.output[0]) + ".chunk.log")
	params:
		chunksize = 1000,
		outdir = os.path.join(TEMP_DIR, "chunked_proteins")
	threads:
		HPC_CONFIG.get_cores("gmc_chunk_proteins")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_chunk_proteins") * attempt
	shell:
		"mkdir -p {params.outdir}" + \
		# awk script by Pierre Lindenbaum https://www.biostars.org/p/13270/
		" && awk 'BEGIN {{n=0;m=1;}} /^>/ {{ if (n%{params.chunksize}==0) {{f=sprintf(\"{params.outdir}/chunk-%d.txt\",m); m++;}}; n++; }} {{ print >> f }}' {input[0]} &> {log}"

rule gmc_metrics_blastp_chunked:
	input:
		chunk = os.path.join(TEMP_DIR, "chunked_proteins", "chunk-{chunk}.txt"),
		db = rules.gmc_metrics_blastp_mkdb.output[0]
	output:
		os.path.join(TEMP_DIR, "chunked_proteins", "{run}", "chunk-{chunk}." + config["blast-mode"] + ".tsv")
	log:
		os.path.join(LOG_DIR, "blast_logs", "chunk-{chunk}.{run}." + config["blast-mode"] + ".log")
	params:
		program_call = config["program_calls"]["blast"].format(program=config["blast-mode"]),
		program_params = config["params"]["blast"][config["blast-mode"]]
	threads:
		HPC_CONFIG.get_cores("gmc_metrics_blastp_chunked")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_metrics_blastp_chunked") * attempt
	shell:
		"{params.program_call} {params.program_params} -query {input.chunk} -out {output[0]} -num_threads {threads} " + \
		"-db {input.db} -outfmt \"6 qseqid sseqid pident qstart qend sstart send qlen slen length nident mismatch positive gapopen gaps evalue bitscore\" &> {log}"


def aggregate_blastp_input(wildcards):
	checkpoint_output = checkpoints.gmc_chunk_proteins.get(**wildcards).output.chunk_dir
	return expand(
		os.path.join(TEMP_DIR, "chunked_proteins", "{run}", "chunk-{chunk}." + config["blast-mode"] + ".tsv"),
		run=wildcards.run,
		chunk=glob_wildcards(os.path.join(checkpoint_output, "chunk-{chunk}.txt")).chunk
	)


rule gmc_metrics_blastp_combine:
	input:
		aggregate_blastp_input
	output:
		os.path.join(EXTERNAL_METRICS_DIR, config["blast-mode"], "{run}", "{run}." + config["blast-mode"] + ".tsv")		
	threads:
		HPC_CONFIG.get_cores("gmc_metrics_blastp_combine")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_metrics_blastp_combine") * attempt
	run:
		with open(output[0], "w") as blast_out:
			for f in input:
				print(open(f).read(), end="", flush=True, file=blast_out)
				os.remove(f)


rule gmc_metrics_blastp_tophit:
	input:
		rules.gmc_metrics_blastp_combine.output[0]
	output:
		rules.gmc_metrics_blastp_combine.output[0] + ".tophit"
	params:
		pident_threshold = config["params"]["blast"]["tophit"]["pident_threshold"],
		qcov_threshold = config["params"]["blast"]["tophit"]["qcov_threshold"]
	threads:
		HPC_CONFIG.get_cores("gmc_metrics_blastp_tophit")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_metrics_blastp_tophit") * attempt
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


rule gmc_metrics_generate_metrics_info:
	input:
		prev_outputs = OUTPUTS,
		bproteins = os.path.join(EXTERNAL_METRICS_DIR, "busco_proteins", "busco_proteins.tsv")
	output:
		os.path.join(EXTERNAL_METRICS_DIR, "metrics_info.txt")
	run:
		import os
		import glob
		import csv

		# this block is unstable against tempering with the structure of the generate_metrics dir
		# might work better as state-machine
		with open(output[0], "wt") as metrics_info:
			walk = os.walk(EXTERNAL_METRICS_DIR, followlinks=True)
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
				elif cwd_base == "repeats":
					mclass = "repeat"
					for path in glob.glob(os.path.join(cwd, "*.no_strand.exon.gff.cbed.parsed.txt")):
						mid = os.path.basename(path).split(".")[0]
						rows.append((ExternalMetrics.REPEAT_ANNOTATION, mclass, mid, path))
				elif cwd_base == "busco_proteins":
					mclass = "busco"
					mid = list(config["data"].get("busco-data", dict()).keys())[0]
					path = os.path.join(EXTERNAL_METRICS_DIR, "busco_proteins", "busco_proteins.tsv")
					rows.append((ExternalMetrics.BUSCO_PROTEINS, mclass, mid, path))
			"""
			for mid in config["data"].get("repeat-data", dict()):
				mclass = "repeat"
				path = config["data"]["repeat-data"][mid][0][0]
				rows.append((ExternalMetrics.REPEAT_ANNOTATION, mclass, mid, path))
			"""

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
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_metrics_generate_metrics_matrix") * attempt
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
		HPC_CONFIG.get_cores("gmc_mikado_serialise")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_mikado_serialise") * attempt
	shell:
		"{params.program_call} {params.program_params} --transcripts {input.transcripts} --external-scores {input.ext_scores} --json-conf {input.config} --procs {threads} -od {params.outdir} &> {log}" + \
		" && touch {output[0]}"

rule gmc_mikado_pick:
	input:
		config = config["mikado-config-file"],
		gtf = rules.gmc_mikado_prepare.output[1],
		serialise_done = rules.gmc_mikado_serialise.output[0],
		db = rules.gmc_mikado_serialise.output[1]
	output:
		loci = os.path.join(config["outdir"], "mikado.loci.gff3"),
		subloci = os.path.join(config["outdir"], "mikado.subloci.gff3")
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="pick"),
		program_params = config["params"]["mikado"]["pick"],
		outdir = config["outdir"]
	threads:
		HPC_CONFIG.get_cores("gmc_mikado_pick")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_mikado_pick") * attempt
	shell:
		"{params.program_call} {params.program_params} -od {params.outdir} --procs {threads} --json-conf {input.config} --subloci-out $(basename {output.subloci}) -db {input.db} {input.gtf}"

rule gmc_parse_mikado_pick:
	input:
		loci = rules.gmc_mikado_pick.output[0]
	output:
		gff = os.path.join(config["outdir"], POST_PICK_PREFIX + ".gff")
	shell:
		"parse_mikado_gff {input.loci} > {output.gff}"

rule gmc_gffread_extract_sequences_post_pick:
	input:
		gff = rules.gmc_parse_mikado_pick.output[0],
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


rule gmc_calculate_cds_lengths_post_pick:
	input:
		rules.gmc_gffread_extract_sequences_post_pick.output.cds
	output:
		rules.gmc_gffread_extract_sequences_post_pick.output.cds + ".lengths"
	params:
		min_cds_length = config["misc"]["min_cds_length"]
	run:
		def read_fasta(f):
			sid, seq = None, list()
			for line in open(f):
				line = line.strip()
				if line.startswith(">"):
					if seq:
						yield sid, "".join(seq)
						seq = list()
					sid = line[1:].split()[0]
				else:
					seq.append(line)
			if seq:
				yield sid, "".join(seq).upper().replace("U", "T")
		with open(output[0], "w") as fout:
			for sid, seq in read_fasta(input[0]):
				has_stop = any(seq.endswith(stop_codon) for stop_codon in ("TAA", "TGA", "TAG"))
				discard = (has_stop and len(seq) < params.min_cds_length) or (not has_stop and len(seq) < params.min_cds_length - 3)
				print(sid, int(discard), sep="\t", file=fout)


rule gmc_gff_genometools_check_post_pick:
	input:
		rules.gmc_parse_mikado_pick.output[0]
	output:
		gff = os.path.join(config["outdir"], POST_PICK_PREFIX + ".gt_checked.gff")
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
		os.path.join(config["outdir"], POST_PICK_PREFIX + ".gt_checked.validation_report.txt")
	threads:
		HPC_CONFIG.get_cores("gmc_gff_validate_post_gt")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_gff_validate_post_gt") * attempt
	shell:
		"validate_gff3 {input} > {output}"

rule gmc_kallisto_index_post_pick:
	input:
		rules.gmc_gffread_extract_sequences_post_pick.output.cdna
	output:
		os.path.join(config["outdir"], "kallisto", os.path.basename(rules.gmc_gffread_extract_sequences_post_pick.output.cdna) + ".idx")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".kallisto_index_post_pick.log")
	params:
		program_call = config["program_calls"]["kallisto"].format(program="index"),
		program_params = config["params"].get("kallisto", {}).get("index", "")
	threads:
		HPC_CONFIG.get_cores("gmc_kallisto_index_post_pick")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_kallisto_index_post_pick") * attempt
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
		HPC_CONFIG.get_cores("gmc_kallisto_quant_post_pick")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_kallisto_quant_post_pick") * attempt
	shell:
		"{params.program_call} {params.program_params} {params.stranded} -i {input.index} -o {params.outdir} --threads {threads} {input.reads} &> {log}"

rule gmc_collapse_metrics:
	input:
		gff = rules.gmc_parse_mikado_pick.output[0],
		ext_scores = rules.gmc_metrics_generate_metrics_matrix.output[0],
		metrics_info = rules.gmc_metrics_generate_metrics_info.output[0],
		expression = expand(rules.gmc_kallisto_quant_post_pick.output, run=config["data"]["expression-runs"].keys()),
		cds_lengths = rules.gmc_calculate_cds_lengths_post_pick.output[0]
	output:
		os.path.join(config["outdir"], POST_PICK_PREFIX + ".collapsed_metrics.tsv"),
		os.path.join(config["outdir"], "COLLAPSE_METRICS_DONE"),
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".collapse_metrics.log")
	threads:
		HPC_CONFIG.get_cores("gmc_collapse_metrics")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_collapse_metrics") * attempt
	run:
		from gmc.scripts.collapse_metrics import MetricCollapser
		mc = MetricCollapser(input.gff, input.metrics_info, input.ext_scores, input.cds_lengths, input.expression)
		with open(output[0], "w") as out:
			mc.write_scores(config["collapse_metrics_thresholds"], stream=out)
		open(output[1], "w").close()

rule gmc_create_release_gffs:
	input:
		gff = rules.gmc_gff_genometools_check_post_pick.output[0],
		metrics_info = rules.gmc_collapse_metrics.output[0],
		sentinel = rules.gmc_collapse_metrics.output[1]
	output:
		os.path.join(config["outdir"], POST_PICK_PREFIX + ".release.unsorted.gff3"),
		os.path.join(config["outdir"], POST_PICK_PREFIX + ".release_browser.unsorted.gff3")
	params:
		annotation_version = config.get("annotation_version", "EIv1"),
		genus_identifier = config.get("genus_identifier", "XYZ")
	threads:
		HPC_CONFIG.get_cores("gmc_create_release_gffs")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_create_release_gffs") * attempt
	shell:
		"create_release_gff3 {input.gff} {input.metrics_info} --annotation-version {params.annotation_version} --genus-identifier {params.genus_identifier} 2> {LOG_DIR}/create_release_gff.log"

rule gmc_sort_release_gffs:
	input:
		rules.gmc_create_release_gffs.output
	output:
		os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".release.gff3"),
		os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".release_browser.gff3")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".sort_release_gffs.log")
	params:
		program_call = config["program_calls"]["genometools"],
		program_params = config["params"]["genometools"]["sort"]
	threads:
		HPC_CONFIG.get_cores("gmc_sort_release_gffs")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_sort_release_gffs") * attempt
	shell:
		"{params.program_call} {params.program_params} {input[0]} > {output[0]} 2> {log}" + \
		" && {params.program_call} {params.program_params} {input[1]} > {output[1]} 2>> {log}"

rule gmc_final_sanity_check:
	input:
		rules.gmc_sort_release_gffs.output[0]
	output:
		os.path.join(RESULTS_DIR, RELEASE_PREFIX + ".sanity_checked.release.gff3")
	log:
		os.path.join(LOG_DIR, config["prefix"] + ".final_sanity_check.log")
	threads:
		HPC_CONFIG.get_cores("gmc_final_sanity_check")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_final_sanity_check") * attempt
	shell:
		"sanity_check {input[0]} > {output[0]} 2> {log}"

rule gmc_mikado_pick_extract_coords:
	input:
		rules.gmc_sort_release_gffs.output[0]
	output:
		rules.gmc_sort_release_gffs.output[0] + ".coords"
	run:
		import csv
		import re
		with open(output[0], "w") as coords_out:
			for row in csv.reader(open(input[0]), delimiter="\t"):
				if row and not row[0].startswith("#") and "rna" in row[2].lower():
					tid = re.search('ID\s?=\s?"?([^";]+)"?;', row[8]).group(1)
					print(tid, "{seq}:{start}..{end}".format(seq=row[0], start=row[3], end=row[4]), sep="\t", flush=True, file=coords_out)

rule gmc_generate_mikado_stats:
	input:
		rules.gmc_final_sanity_check.output[0]
	output:
		rules.gmc_final_sanity_check.output[0] + ".mikado_stats.txt",
		rules.gmc_final_sanity_check.output[0] + ".mikado_stats.tsv"
	params:
		program_call = config["program_calls"]["mikado"].format(container=config["mikado-container"], program="util stats"),
	threads:
		HPC_CONFIG.get_cores("gmc_generate_mikado_stats")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_generate_mikado_stats") * attempt
	shell:
		"{params.program_call} {input} --tab-stats {output[1]} > {output[0]}" + \
		" && parse_mikado_stats {output[0]} > {output[0]}.summary"

rule gmc_extract_final_sequences:
	input:
		gff = rules.gmc_final_sanity_check.output[0],
		refseq = config["reference-sequence"]
	output:
		cdna = rules.gmc_final_sanity_check.output[0] + ".cdna.fasta",
		tbl = rules.gmc_final_sanity_check.output[0] + ".gffread.table.txt",
		cds = rules.gmc_final_sanity_check.output[0] + ".cds.fasta",
		pep = rules.gmc_final_sanity_check.output[0] + ".pep.raw.fasta"
	params:
		program_call = config["program_calls"]["gffread"],
		table_format = "--table @chr,@start,@end,@strand,@numexons,@covlen,@cdslen,ID,Note,confidence,representative,biotype,InFrameStop,partialness"
	shell:
		"{params.program_call} {input.gff} -g {input.refseq} -P {params.table_format} -W -w {output.cdna} -x {output.cds} -y {output.pep} -o {output.tbl}"

rule gmc_cleanup_final_proteins:
	input:
		rules.gmc_extract_final_sequences.output.pep
	output:
		rules.gmc_extract_final_sequences.output.pep.replace(".raw.fasta", ".fasta")
	log:
		os.path.join(LOG_DIR, "cleanup_proteins.log")
	params:
		prefix = rules.gmc_extract_final_sequences.output.pep.replace(".raw.fasta", ""),
		program_call = config["program_calls"]["prinseq"],
		program_params = config["params"]["prinseq"]
	threads:
		HPC_CONFIG.get_cores("gmc_cleanup_final_proteins")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_cleanup_final_proteins") * attempt
	shell:
		"{params.program_call} -aa -fasta {input} {params.program_params} -out_good {params.prefix} -out_bad {params.prefix}.bad"

rule gmc_generate_full_table:
	input:
		stats_table = rules.gmc_generate_mikado_stats.output[1],
		seq_table = rules.gmc_extract_final_sequences.output.tbl,
		bt_conf_table = rules.gmc_final_sanity_check.output[0]
	output:
		final_table = rules.gmc_final_sanity_check.output[0] + ".final_table.tsv",
		summary = rules.gmc_final_sanity_check.output[0] + ".biotype_conf.summary"
	threads:
		HPC_CONFIG.get_cores("gmc_generate_full_table")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("gmc_generate_full_table") * attempt
	run:
		import csv
		from collections import Counter
		colheaders = ["Confidence", "Biotype", "InFrameStop", "Partialness"]
		pt_cats = {".": "complete", "5_3": "fragment", "3": "3prime_partial", "5": "5prime_partial"}
		if_pt = dict((row[7], row[12:14]) for row in csv.reader(open(input.seq_table), delimiter="\t"))
		genes, transcripts = dict(), dict()
		for row in csv.reader(open(input.bt_conf_table), delimiter="\t"):
			if not row[0].startswith("#"):
				if row[2].lower() in {"gene", "mrna", "ncrna"}:
					attrib = dict(item.split("=") for item in row[8].strip().split(";"))
					if row[2].lower() == "gene":
						genes[attrib["ID"]] = (attrib["biotype"], attrib["confidence"])
					else:
						transcripts[attrib["ID"]] = (genes.get(attrib["Parent"], (None, None))[0], attrib["confidence"])

		r = csv.reader(open(input.stats_table), delimiter="\t")
		head = ["#{}.{}".format(c, col) for c, col in enumerate(next(r), start=1)]
		head.extend("#{}.{}".format(c, col) for c, col in enumerate(colheaders, start=len(head)+1))
		with open(output.final_table, "w") as tbl_out:
			print(*head, sep="\t", flush=True, file=tbl_out)
			for row in r:
				row.extend(transcripts.get(row[0], [".", "."])[::-1])
				row.extend(if_pt.get(row[0], [".", "."]))
				row[-1] = pt_cats.get(row[-1], "NA")
				print(*row, sep="\t", flush=True, file=tbl_out)

		tcounts, gcounts = Counter(transcripts.values()), Counter(genes.values())
		with open(output.summary, "w") as summary_out:
			print("Biotype", "Confidence", "Gene", "Transcript", sep="\t", file=summary_out)
			categories = set(tcounts).union(set(gcounts))
			for cat in sorted(categories, key=lambda x:gcounts[x], reverse=True):
				print(*cat, gcounts[cat], tcounts[cat], sep="\t", file=summary_out)
			print("Total", "", sum(gcounts.values()), sum(tcounts.values()), sep="\t", file=summary_out)

rule split_proteins_prepare:
	input:
		rules.gmc_gffread_extract_sequences.output[0]
	output:
		expand(os.path.join(BUSCO_PATH, "runs", "proteins_prepare", "input", "{run}.proteins.fasta"), run=config["transcript_models"])
	log:
		os.path.join(BUSCO_PATH, "logs", "split_proteins_prepare.log")
	run:
		from gmc.scripts.busco_splitter import split_fasta
		fasta_files = {tm: open(os.path.join(BUSCO_PATH, "runs", "proteins_prepare", "input", tm + ".proteins.fasta"), "w") for tm in config["transcript_models"]}
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
		copy = config["program_calls"]["copy"]
	threads:
		HPC_CONFIG.get_cores("busco_proteins_prepare")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("busco_proteins_prepare") * attempt
	shell:
		BUSCO_CMD

rule busco_proteins_final:
	input:
		rules.gmc_extract_final_sequences.output.pep
	output:
		os.path.join(BUSCO_PATH, "runs", "proteins_final", "proteins_final", "run_{}".format(BUSCO_LINEAGE), "short_summary.txt"),
		os.path.join(BUSCO_PATH, "runs", "proteins_final", "proteins_final", "run_{}".format(BUSCO_LINEAGE), "full_table.tsv"),
		os.path.join(BUSCO_PATH, "runs", "proteins_final", "proteins_final", "run_{}".format(BUSCO_LINEAGE), "missing_busco_list.tsv")
	log:
		os.path.join(BUSCO_PATH, "logs", "proteins_final.log")
	params:
		input = os.path.abspath(rules.gmc_extract_final_sequences.output.pep),
		program_call = config["program_calls"]["busco"],
		program_params = config["params"]["busco"]["proteins_final"],
		lineage_path = config["busco_analyses"]["lineage"],
		run = "proteins_final",
		busco_mode = "proteins",
		busco_stage = "proteins_final",
		copy = config["program_calls"]["copy"]
	threads:
		HPC_CONFIG.get_cores("busco_proteins_final")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("busco_proteins_final") * attempt
	shell:
		BUSCO_CMD


rule split_transcripts_prepare:
	input:
		rules.gmc_gffread_extract_sequences.output[1]
	output:
		expand(os.path.join(BUSCO_PATH, "runs", "transcripts_prepare", "input", "{run}.cdna.fasta"), run=config["transcript_models"])
	log:
		os.path.join(BUSCO_PATH, "logs", "split_transcripts_prepare.log")
	run:
		from gmc.scripts.busco_splitter import split_fasta
		fasta_files = {tm: open(os.path.join(BUSCO_PATH, "runs", "transcripts_prepare", "input", tm + ".cdna.fasta"), "w") for tm in config["transcript_models"]}
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
		copy = config["program_calls"]["copy"]
	threads:
		HPC_CONFIG.get_cores("busco_transcripts_prepare")
	resources:
		mem_mb = lambda wildcards, attempt: HPC_CONFIG.get_memory("busco_transcripts_prepare") * attempt
	shell:
		BUSCO_CMD

rule busco_transcripts_final:
	input:
		rules.gmc_extract_final_sequences.output.cdna
	output:
		os.path.join(BUSCO_PATH, "runs", "transcripts_final", "transcripts_final", "run_{}".format(BUSCO_LINEAGE), "short_summary.txt"),
		os.path.join(BUSCO_PATH, "runs", "transcripts_final", "transcripts_final", "run_{}".format(BUSCO_LINEAGE), "full_table.tsv"),
		os.path.join(BUSCO_PATH, "runs", "transcripts_final", "transcripts_final", "run_{}".format(BUSCO_LINEAGE), "missing_busco_list.tsv"),
	log:
		os.path.join(BUSCO_PATH, "logs", "transcripts_final.log")
	params:
		input = os.path.abspath(rules.gmc_extract_final_sequences.output.cdna),
		program_call = config["program_calls"]["busco"],
		program_params = config["params"]["busco"]["transcripts_final"],
		lineage_path = config["busco_analyses"]["lineage"],
		run = "transcripts_final",
		busco_mode = "transcriptome",
		busco_stage = "transcripts_final",
		copy = config["program_calls"]["copy"]
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
			copy = config["program_calls"]["copy"]
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
		rules.gmc_mikado_prepare_extract_coords.output[0],
		rules.gmc_mikado_pick_extract_coords.output[0],
		BUSCO_ANALYSES + BUSCO_PROTEIN_PREPARE_RUNS
	output:
		BUSCO_TABLE
	run:
		from gmc.scripts.generate_busco_tables import BuscoTableGenerator

		btg = BuscoTableGenerator(
			os.path.join(config["outdir"], "tx2gene"),
			input[0],
			input[1],
			os.path.join(BUSCO_PATH, "runs")	
		)
		btg.write_review_table(output[0])
		btg.write_raw_data(output[0])
		btg.write_busco_table(output[0], config["misc"]["busco_max_copy_number"])	
