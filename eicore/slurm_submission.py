#!/usr/bin/env python

import os
import sys
import argparse
import yaml
import subprocess
from abc import ABC, abstractmethod

from eicore.snakemake_helper import DEFAULT_CONFIG_FILE, NOW

DEFAULT_PARTITION = "tgac-pap"
DEFAULT_EMAIL = "$(whoami)@nbi.ac.uk"

SBATCH_TEMPLATE = "-J {program}_{ticket}_{dataset_id} -o {dataset_id}.{short_program}{resume_str}.out.log -e {dataset_id}.{short_program}{resume_str}.err.log --mail-type=FAIL --mail-user={email} -p {partition}"

# https://stackoverflow.com/questions/431684/how-do-i-change-directory-cd-in-python
from contextlib import contextmanager


@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


class SlurmSubmitter(ABC):
    @staticmethod
    def generate_resume_string(output_dir, prefix):
        output_dir = os.path.basename(output_dir)

        def resume_files_exist(prefix, r):
            return any(
                os.path.exists(
                    "{prefix}.resume{attempt}.{logtype}.log".format(
                        prefix=prefix, attempt=r, logtype=logtype
                    )
                )
                for logtype in ("out", "err")
            )

        file_prefix = "{output_dir}.{prefix}".format(
            output_dir=output_dir, prefix=prefix
        )
        r = 1
        while resume_files_exist(file_prefix, r):
            r += 1

        return "{}.{}".format(file_prefix, r)

    @abstractmethod
    def compile_slurm_options(self):
        pass

    def __init__(self, cmd_args=sys.argv[1:]):
        print("in SlurmSubmitter")
        print(cmd_args)
        if not hasattr(self, "whoami"):
            raise ValueError("Error: Cannot run without specifying module.")

        config = yaml.load(open(DEFAULT_CONFIG_FILE), Loader=yaml.SafeLoader)
        self.argparser = self.argparser_func(config)
        self.argparser.add_argument("--no-submit", action="store_true")
        args = self.argparser.parse_args(cmd_args)

        self.email = os.environ.get("SNAKEMAKE_EMAIL", DEFAULT_EMAIL)
        self.partition = (
            args.partition
            if hasattr(args, "partition") and args.partition
            else DEFAULT_PARTITION
        )

        self.__compile_program_options(
            args, ignore_program_options=self.ignore_program_options
        )
        self.require_resume_mode = True

    def generate_submission_string(self):
        #  adding output_dir -o, as it was ignored in eipap_option parsing: easier to deal with custom output_dir
        od_string = (" -o {}".format(self.output_dir)) if self.output_dir else ""
        self.submission_string = "{whoami}{output_dir} {program_options}".format(
            whoami=self.whoami,
            output_dir=od_string,
            program_options=" ".join(self.program_options),
        )

    def __remove_submission_script(self):
        try:
            os.remove(self.submission_script)
        except:
            pass

    def __write_submission_script(self):
        try:
            with open(self.submission_script, "w") as cmd_out:
                print("#!/bin/bash -e", self.submission_string, sep="\n", file=cmd_out)
        except:
            raise ValueError("Could not write submission script.")

    def __compile_program_options(self, args, ignore_program_options):
        self.program_options = list()
        for k, v in vars(args).items():
            if k not in ignore_program_options and v is not None and v:
                opt = "--" + k
                if not type(v) is bool:
                    opt += "=" + str(v)
                self.program_options.append(opt)

        try:
            self.jira = args.jira
        except:
            self.jira = None

        try:
            self.output_dir = args.output_dir
        except:
            self.output_dir = ""

        # for some reason this was commented out (one of the pap tools?)
        if self.output_dir:
            self.program_options.append("-o {}".format(self.output_dir))

        self.no_submit = args.no_submit
        try:
            self.program_options.append(args.input)
        except:
            pass
        print(*self.program_options)

    def set_resume_string(self, prefix, program):
        if not self.require_resume_mode or (
            hasattr(self, "mode") and self.mode.upper() == "RESUME"
        ):
            self.resume_string = SlurmSubmitter.generate_resume_string(prefix, program)
        else:
            self.resume_string = ""

    def run(self, cleanup=True):

        with cd(self.rundir if self.rundir is not None else "."):
            if self.jira is not None:
                print("Changing into QC directory for {0} ".format(self.jira))
            else:
                print("Changing into rundir {} ".format(self.rundir))

            self.__remove_submission_script()
            self.__write_submission_script()

            print(
                'Running submission command: sbatch {slurm_opts} --wrap="{whoami} {prog_opts}"'.format(
                    whoami=self.whoami,
                    slurm_opts=self.slurm_options,
                    prog_opts=" ".join(map(str, self.program_options)),
                )
            )

            if not self.no_submit:
                p = subprocess.Popen(
                    "sbatch {} {}".format(self.slurm_options, self.submission_script),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    shell=True,
                )

                o, e = p.communicate()

                print(o.decode())
                if e.decode():
                    print(e.decode())

            if cleanup:
                self.__remove_submission_script()


class PapToolSubmitter(SlurmSubmitter):
    @abstractmethod
    def compile_slurm_options(self):
        pass

    def determine_run_directory(self):
        print("Determining run directory" + self.jira)
        try:
            p = subprocess.Popen(
                "eipap_cd qc {}".format(self.jira),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.rundir = p.communicate()[0].decode().strip().split("\n")[-1]
        except:
            raise ValueError("Cannot determine qc directory for {}".format(self.jira))

    def __init__(self, cmd_args=sys.argv[1:]):
        print("in PapToolSubmitter")
        self.default_output_dir = "PAP_" + NOW.split("_")[0]
        self.ignore_program_options = {"input", "output_dir", "no_submit"}
        super().__init__(cmd_args=cmd_args)

