# coding: utf-8
from setuptools import setup, find_packages
from setuptools.extension import Extension
from distutils.extension import Extension
from codecs import open
from os import path
import glob
import re
import sys


here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
	description = long_description = description.read()

name="gmc"
version = "1.2.2"

if sys.version_info.major != 3:
	raise EnvironmentError("""gmc is a python module that requires python3, and is not compatible with python2. Also, it is now 2020 and support for 2.x has ceased.""")


setup(
	name=name,
	version=version,
	description=description,
	long_description=long_description,
	author="Christian Schudoma",
	author_email="christian.schudoma@earlham.ac.uk",
	license="MIT",
	classifiers=[
		"Development Status :: 4 - Beta",
		"Topic :: Scientific Engineering :: Bio/Informatics",
		"License :: OSI Approved :: MIT License",
		"Operating System :: POSIX :: Linux",
		'Programming Language :: Python :: 3.4',
		"Programming Language :: Python :: 3.5",
		"Programming Language :: Python :: 3.6"
		"Programming Language :: Python :: 3.7"
	],
	zip_safe=False,
	keywords="gene annotation",
	packages=find_packages(exclude=["test"]),
	scripts=[
		script for script in glob.glob("bin/slurm/*_sub")
	],
	install_requires=[
		"snakemake>5.4.0",
		"drmaa",
	],
	entry_points={
		"console_scripts": [
			"gmc=gmc.__main__:main",
			"generate_metrics=gmc.scripts.generate_metrics:main",
			"parse_mikado_gff=gmc.scripts.parse_mikado_gff:main",
			"protein_completeness=gmc.scripts.protein_completeness:main",
			"collapse_metrics=gmc.scripts.collapse_metrics:main",
			"validate_gff3=gmc.scripts.validate_gff3:main",
			"create_release_gff3=gmc.scripts.create_release_gff:main",
			"sanity_check=gmc.scripts.sanity_check:main",
			"parse_mikado_stats=gmc.scripts.parse_mikado_stats:main",
			"analyse_busco=gmc.scripts.analyse_busco:main"
		]
	},
	package_data={
		"gmc.zzz": ["*.smk"]
	},
	include_package_data=True,
	data_files=[
		("etc", glob.glob("etc/*.*"))
	]
)

