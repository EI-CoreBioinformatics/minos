import os

localrules: all

rule all:
	input:
		[os.path.join(config["outdir"], "mikado_prepared.fasta"), os.path.join(config["outdir"], "generate_metrics", "CPC-2.0_beta", "mikado_prepared.fasta.cpc2output.txt")]

rule gmc_mikado_prepare:
	input:
		config["mikado-config-file"]
	output:
		# os.path.join(config["outdir"], "MIKADO_PREPARE_DONE")
		os.path.join(config["outdir"], "mikado_prepared.fasta"),
		os.path.join(config["outdir"], "mikado_prepared.gtf")
	params:
		mikado = config["mikado-container"] + " mikado prepare",
		min_length = 100,
		outdir = config["outdir"]
	log:
		os.path.join(config["outdir"], "logs", config["prefix"] + ".mikado_prepare.log")
	threads:
		30
	shell:
		"singularity exec {params.mikado} --minimum_length {params.min_length} --json-conf {input[0]} --procs {threads} -od {params.outdir} &> {log}" #fasta && touch {output[0]}"


rule gmc_metrics_cpc2:
	input:
		rules.gmc_mikado_prepare.output[0]
	output:
		#rules.gmc_mikado_prepare.output[0] + ".cpc2output.txt"
		os.path.join(config["outdir"], "generate_metrics", "CPC-2.0_beta", os.path.basename(rules.gmc_mikado_prepare.output[0]) + ".cpc2output.txt")
	params:
		x = 1
	log:
		os.path.join(config["outdir"], "logs", config["prefix"] + ".CPC2.log")
	threads:
		4
	shell:
		"set +u && source CPC-2.0_beta_py3_cs && " + \
		"/usr/bin/time -v CPC2.py -r -i {input[0]} -o {output[0]} &> {log} "
			



#		sbatch -p ei-medium -c 4 --mem 20G -o out_CPC2.py.%j.log -J CPC --wrap "ln -s /ei/workarea/group-ga/Projects/CB-GENANNO-430_Annotation_of_Melia_azedarach_and_Quillaja_saponaria/Analysis/Quillaja_saponaria/mikado-1.2.4/annotation_run1/run1/mikado_prepared.fasta && source CPC-2.0_beta && /usr/bin/time -v CPC2.py -r -i mikado_prepared.fasta -o mikado_prepared.fasta.cpc2output.txt"

# /ei/workarea/group-ga/Projects/CB-GENANNO-430_Annotation_of_Melia_azedarach_and_Quillaja_saponaria/Analysis/Quillaja_saponaria/mikado-1.2.4/annotation_run1/run1/generate_metrics/CPC-2.0_beta/mikado_prepared.fasta.cpc2output.txt

