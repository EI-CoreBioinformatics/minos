import os

localrules: all, gmc_metrics_leaff

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

OUTPUTS = [
	os.path.join(config["outdir"], "mikado_prepared.fasta"), 
	os.path.join(EXTERNAL_METRICS_DIR, "CPC-2.0_beta", "mikado_prepared.fasta.cpc2output.txt"),
	os.path.join(config["outdir"], "logs", "mikado_prepared.fasta.leaff.log"),
	os.path.join(EXTERNAL_METRICS_DIR, "kallisto", "mikado_prepared.fasta.idx")
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


rule all:
	input:
		OUTPUTS

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

rule gmc_metrics_leaff:
	input:
		rules.gmc_mikado_prepare.output[0]
	output:
		os.path.join(LOGDIR, os.path.basename(rules.gmc_mikado_prepare.output[0]) + ".leaff.log")
	params:
		chunksize = 500,
		outdir = os.path.join(config["outdir"], "tmp"),
		wrapper = "/ei/workarea/group-ga/Scripts/leaff_v0.2.pl"
	shell:
		"set +u && source perl-5.20.1_gk && mkdir -p {params.outdir} && " + \
		"{params.wrapper} {input[0]} {params.outdir}/chunk {params.chunksize} &> {output[0]}"

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
		"singularity exec {params.mikado} --extended-refmap -r {input.mika} -p {input.transcripts} -o vs_transcripts_{params.transcripts} &> {log} && " + \
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
		"singularity exec {params.mikado} --exclude-utr --extended-refmap -r {input.mika} -p {input.proteins} -o vs_proteins_{params.proteins} &> {log} && " + \
		"touch {output[0]}"

#rule gmc_metrics_blastp:
