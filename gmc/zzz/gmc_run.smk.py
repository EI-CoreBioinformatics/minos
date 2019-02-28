import os

localrules: all

rule all:
	input:
		[os.path.join(config["outdir"], "MIKADO_PREPARE_DONE")]

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



#WD: /ei/workarea/group-ga/Projects/CB-GENANNO-430_Annotation_of_Melia_azedarach_and_Quillaja_saponaria/Analysis/Quillaja_saponaria/mikado-1.2.4/annotation_run1/run1
#CMD:
#sbatch -p ei-medium -c 30 --mem 120G -o out_mikado.prepare.%j.log -J QS_Mikado_prepare --wrap "source mikado-devel && /usr/bin/time -v mikado prepare --minimum_length 100 --procs 30 --json-conf Quillaja_saponaria.configuration.yaml"

		
