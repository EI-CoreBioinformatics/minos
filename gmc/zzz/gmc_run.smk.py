import os

localrules: all, gmc_metrics_leaff

EXTERNAL_METRICS_DIR = os.path.join(config["outdir"], "generate_metrics")
LOGDIR = os.path.join(config["outdir"], "logs")




def get_rnaseq(wc):
	# return " ".join(" ".join(pair) for pair in config["expression-runs"][wc.run])
	return [item for sublist in config["expression-runs"][wc.run] for item in sublist]


OUTPUTS = [
	os.path.join(config["outdir"], "mikado_prepared.fasta"), 
	os.path.join(EXTERNAL_METRICS_DIR, "CPC-2.0_beta", "mikado_prepared.fasta.cpc2output.txt"),
	os.path.join(config["outdir"], "logs", "mikado_prepared.fasta.leaff.log"),
	os.path.join(EXTERNAL_METRICS_DIR, "kallisto", "mikado_prepared.fasta.idx")
]

for run in config.get("expression-runs", {}):
	OUTPUTS.append(os.path.join(EXTERNAL_METRICS_DIR, "kallisto", run, "abundance.tsv"))



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
		stranded = "--rf-stranded" if True else "", # TODO!
		outdir = os.path.join(EXTERNAL_METRICS_DIR, "kallisto", "{run}")
	threads:
		32
	shell:
		"set +u && source kallisto-0.44.0 && /usr/bin/time -v kallisto quant {params.stranded} -i {input.index} -o {params.outdir} -b 100 --threads {threads} {input.reads} &> {log}"



"""
id: SRA_kallisto
      type: expression
      multiplier: 0
      not_fragmentary_min_value: 0
      files: [[/ei/workarea/group-pb/CB-PPBFX-550_Anne_Osbourn_JIC_AO_ENQ-2696_A_01_ext/Reads/ERR706840_1.fastq.gz, /ei/workarea/group-pb/CB-PPBFX-550_Anne_Osbourn_JIC_AO_ENQ-2696_A_01_ext/Reads/ERR706840_2.fastq.gz]]
    -
"""                                                                                                                                                                                                                      			

#		sbatch -p ei-medium -c 4 --mem 20G -o out_CPC2.py.%j.log -J CPC --wrap "ln -s /ei/workarea/group-ga/Projects/CB-GENANNO-430_Annotation_of_Melia_azedarach_and_Quillaja_saponaria/Analysis/Quillaja_saponaria/mikado-1.2.4/annotation_run1/run1/mikado_prepared.fasta && source CPC-2.0_beta && /usr/bin/time -v CPC2.py -r -i mikado_prepared.fasta -o mikado_prepared.fasta.cpc2output.txt"

# /ei/workarea/group-ga/Projects/CB-GENANNO-430_Annotation_of_Melia_azedarach_and_Quillaja_saponaria/Analysis/Quillaja_saponaria/mikado-1.2.4/annotation_run1/run1/generate_metrics/CPC-2.0_beta/mikado_prepared.fasta.cpc2output.txt

