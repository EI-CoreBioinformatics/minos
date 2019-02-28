import os

localrules: all, gmc_mikado_prepare

rule all:
	input:
		[os.path.join(config["outdir"], "DUMMY")]

rule gmc_mikado_prepare:
	input:
		config["mikado-config-file"]
	output:
		os.path.join(config["outdir"], "DUMMY")
	params:
		mikado = config["mikado-container"] + " mikado prepare"
	shell:
		"singularity exec {params.mikado} -h && touch {output[0]}"
		
