import os
import sys
import pkg_resources
import pathlib

LOG_DIR = os.path.join(config["output"]["outdir"], "logs")
ENV_DIR = pkg_resources.resource_filename("minos.zzz", "envs")


localrules: minos_mikado_configure, all

rule all:
	input:
		config["output"]["mikado-config-file"]

# MIKADO_CONFIGURE_CMD = "{cmd} --list {list_file}{external_metrics}-od {output_dir} --reference {reference} --scoring {scoring_file}{junctions}{mikado_config_file} --full"

rule minos_mikado_configure:
	input:
		tr_list = config["input"]["transcript_models"],
		reference = config["input"]["reference"],
		scoring_template = config["input"]["scoring_template"]
	output:
		config["output"]["mikado-config-file"]
	params:
		program_call = config["program_calls"]["mikado_configure"], 
		outdir = config["output"]["outdir"],
		external_metrics = " --external {} ".format(config["input"]["external_metrics"]) if config["input"]["external_metrics"] else " ",
		junctions = " --junctions {} ".format(config["input"]["junctions"]) if config["input"]["junctions"] else " "
	log:
		os.path.join(LOG_DIR, config["output"]["prefix"] + ".mikado_configure.log")
	conda:
		os.path.join(ENV_DIR, "mikado.yaml")
	shell:
		"{params.program_call} --full" + \
		" --list {input.tr_list}{params.external_metrics}-od {params.outdir} --reference {input.reference}" + \
		" --scoring {input.scoring_template}{params.junctions}{output[0]} &> {log}"

