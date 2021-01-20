import datetime
import os
import time
import sys
from enum import Enum, unique
from textwrap import dedent

from snakemake import snakemake

from .capturing import Capturing


NOW = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H%M%S')

def get_etc_dir():
	path, prefix = __file__, None
	while path:
		path, here = os.path.split(path)
		if here == "lib":
			prefix = path
			break
	
	if prefix is None:
		raise ValueError("Cannot deduce install location from " + __file__)

	return os.path.join(prefix, "../minos/etc")
	

ETC_DIR = get_etc_dir() 
DEFAULT_HPC_CONFIG_FILE = os.path.join(ETC_DIR, "hpc_config.json")
DEFAULT_CONFIG_FILE = os.path.join(ETC_DIR, "minos_config.yaml")

@unique
class RunMode(Enum):
	NORMAL = 0
	RESTART = 1
	RESUME = 2
	VALIDATE = 3
	DRYRUN = 4


def make_exeenv_arg_group(parser, default_hpc_config_file=DEFAULT_HPC_CONFIG_FILE, allow_mode_selection=True, silent=False):  
	"""
	This adds a command line option group to the provided parser.  These options help control the external_process process over various architectures and schedulers
	:param parser: The command line parser to add exeenv options to
	:return: The altered parse
	"""
	if not silent:
		print("Making exe_env arg group")
	hpc_group = parser.add_argument_group("HPC Options",
										  "Controls for how jobs should behave across the HPC resources.")
	if allow_mode_selection:
		hpc_group.add_argument("--mode", choices=[rm.name.lower() for rm in RunMode], default="normal",
							help=dedent("""This option controls how to run the pipelines.
								\"normal\" - This mode will start the run from scratch into an non-existant output directory.  This 
								is the default option.  However, if the output directory already exists then the pipeline will fail.
								\"resume\" - This skips the validation step and resumes the pipeline from an existing output directory 
								where it left off before.
								\"restart\" - Wipes the output directory and start afresh.  WARNING: this will delete all contents in
								the output directory.
								\"validate\" - Does not attempt to run an snakemake pipeline.  Simply analyses and validates any input provided.
								\"dryrun\" - As \'validate\' mode but also performs a dryrun of the pipeline, ensuring the pipeline
								 works with the input as well.  In addition, this produces a DAG in dot format, that can be turned into a plot.  
								 e.g. dot -T svg pipeline.dot > pipeline.svg."""))

	hpc_group.add_argument("--partition", type=str,
						   help="Will run all child jobs on this partition/queue, this setting overrides anything specified in the \"--hpc_config\" file.")
	hpc_group.add_argument("--scheduler", type=str,
						   help="The job scheduler to use.  LSF, PBS and SLURM currently supported.  If running without a scheduler type NONE here. Assumes SLURM by default.")
	hpc_group.add_argument("--no_drmaa", action='store_true', default=False,
						   help="Use this flag if DRMAA is not available")
	hpc_group.add_argument("-N", "--max_nodes", type=int, default=100,
						   help="Maximum number of nodes to use concurrently")
	hpc_group.add_argument("-c", "--max_cores", type=int, default=200,
						   help="Maximum number of cores to use concurrently")
	hpc_group.add_argument("--hpc_config", default=default_hpc_config_file,
						   help="Configuration file for the HPC.  Can be used to override what resources and partitions each job uses.")
	hpc_group.add_argument("--unlock", action='store_true', default=False,
						help="""If the snakemake pipeline is not running because it is reporting that the directory is locked, 
						then you can unlock it using this option.  Please make sure that there are no other snakemake jobs 
						are running in this directory before using this option!""")

	return parser


class ExecutionEnvironment:
	"""Provides some extra functionality for managing various scheduled environments"""
	def __init__(self, args=None, now="eicore", job_suffix=None, log_dir="logs"):
		self.use_drmaa = False
		self.use_scheduler = False
		self.partition = ""
		self.max_nodes = 1
		self.max_cores = 4
		self.sub_cmd = ""
		self.res_cmd = ""
		self.hpc_config = ""
		if args:
			scheduler = args.scheduler if args.scheduler and args.scheduler != '' else "SLURM"
			self.partition = args.partition if args.partition else "{cluster.partition}"
			self.hpc_config = args.hpc_config
			self.use_drmaa = not args.no_drmaa
			self.max_nodes = args.max_nodes
			self.max_cores = args.max_cores
			log_prefix = os.path.join(log_dir, now + "{rule}%j_%N")
			job_name = "{rule}_" + job_suffix if job_suffix and job_suffix != "" else "{rule}"
			if scheduler.upper() == "LSF":
				self.sub_cmd = "bsub"
				self.res_cmd = " -R rusage[mem={cluster.memory}]span[ptile={threads}] -n {threads} -q " + self.partition + " -J " + job_name + " -oo " + log_prefix + ".lsf.log"
				self.use_scheduler = True
			elif scheduler.upper() == "PBS":
				self.sub_cmd = "qsub"
				self.res_cmd = " -lselect=1:mem={cluster.memory}MB:ncpus={threads} -q " + self.partition + " -N " + job_name + " -o " + log_prefix + ".pbs.stdout -e " + log_prefix + ".pbs.stderr"
				self.use_scheduler = True
			elif scheduler.upper() == "SLURM":
				self.sub_cmd = "sbatch"
				self.res_cmd = " -c {cores} -p {partition} --exclude={exclude_nodes} --mem={memory} -J {job_name} -o {logfile} --time={time_limit}".format(
					cores="{threads}",
					partition=self.partition,
					exclude_nodes="{cluster.exclude}",
					memory="{resources.mem_mb}",
					job_name=job_name,
					logfile=log_prefix + ".slurm.log",
					time_limit="{cluster.time}"
				)
				#self.res_cmd =	 " -c {threads}" + \
				#				" -p " + self.partition + \
				#				" --exclude={cluster.exclude}" + \
				#				" --mem={cluster.memory}" + \
				#				" -J " + job_name + \
				#				" -o " + log_prefix + ".slurm.log" + \
				#				" --time={cluster.time}"
				self.use_scheduler = True
			elif scheduler == "" or scheduler.upper() == "NONE":
				pass
			else:
				raise ValueError("Unexpected scheduler configuration.  Check settings.")

	def __str__(self):
		return '\n'.join([
			"Use scheduler: " + str(self.use_scheduler),
			"Use DRMAA: " + str(self.use_drmaa),
			"Submission command: " + self.sub_cmd,
			"Resource command: " + self.res_cmd,
			"Partition: " + self.partition,
			"Max nodes: " + str(self.max_nodes),
			"Max cores: " + str(self.max_cores),
			"HPC configuration file: " + self.hpc_config
		])


def loadPreCmd(*args):
	
	"""Used to prefix a shell command that utilises some external software with another command used
	to load that software"""


	final_command = ""

	for command in args:
		if command:
			cc = command.strip()
			if cc:
				final_command += "set +u && {} &&".format(cc)
			continue

	return final_command


def run_snakemake(snakefile, out_dir, cfg_file, exe_env, dryrun=False, unlock=False):
	"""Helper function for calling external_process pipeline.  This helps us deal with all the different options that we
	might want to switch between running in dryrun and regular mode."""
	res = False
	if dryrun:
		print("Dry run requested.  Will not execute tasks.")
		print()
		with Capturing() as output:
			res = snakemake(snakefile,
							cores=exe_env.max_cores,
							nodes=exe_env.max_nodes,
							configfiles=[cfg_file],
							workdir=".",
							unlock=unlock,
							# force_incomplete=args.force_incomplete,
							# detailed_summary=args.detailed_summary,
							# list_resources=args.list_resources,
							latency_wait=60 if exe_env.use_scheduler else 1,
							printdag=dryrun,
							dryrun=False,
							forceall=dryrun,
							# allowed_rules=args.allowed_rules
							)
		if res:
			dag_file = os.path.join(out_dir,
									"eicore-" + os.path.basename(snakefile).split('.')[0] + "_dag-" + NOW + ".dot")
			print("Saving DAG of pipeline to " + dag_file + " ... ", end="", flush=True)
			with open(dag_file, 'w') as df:
				print("\n".join(output), file=df)
			print("done.")
		else:
			print("Error occured processing Illumina PAP:\n" + "\n".join(output))

	else:

		cluster_cfg = exe_env.hpc_config if exe_env.use_scheduler else None
		cluster = (
			exe_env.sub_cmd + exe_env.res_cmd if not exe_env.use_drmaa else None) if exe_env.use_scheduler else None
		drmaa = exe_env.res_cmd if exe_env.use_drmaa else None
		res = snakemake(snakefile,
						cores=exe_env.max_cores,
						local_cores=exe_env.max_cores,
						nodes=exe_env.max_nodes,
						configfiles=[cfg_file],
						workdir=".",
						cluster_config=cluster_cfg,
						cluster=cluster,
						drmaa=drmaa,
						unlock=unlock,
						printshellcmds=True,
						printreason=True,
						stats=os.path.join(out_dir, os.path.basename(snakefile) + "-" + NOW + ".stats"),
						jobname="minos.{rulename}.{jobid}",
						force_incomplete=True,
						# detailed_summary=args.detailed_summary,
						# list_resources=True,
						latency_wait=60 if exe_env.use_scheduler or exe_env.use_drmaa else 1,
						printdag=dryrun,
						dryrun=dryrun,
						forceall=dryrun,
						verbose=True,
						use_conda=True,
						use_singularity=True,
						restart_times=3,
						max_status_checks_per_second=30,
						keepgoing=True
						# allowed_rules=args.allowed_rules
						)
	return res
