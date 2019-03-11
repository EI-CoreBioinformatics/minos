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
version = "0.1"

if sys.version_info.major != 3:
	raise EnvironmentError("""gmc is a python module that requires python3, and is not compatible with python2. Also, it is now 2019 and support for 2.x will cease soon.""")


setup(
	name=name,
	version=version,
	description=description,
	long_description=long_description,
	# url="https://github.com/EI-CoreBioinformatics/bg",
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
	],
	zip_safe=False,
	keywords="gene annotation",
	packages=find_packages(exclude=["test"]),
	scripts=[
		script for script in glob.glob("gmc/bin/slurm/*_sub")
	],
	install_requires=[
		"snakemake>=5.4.0",
		"drmaa"
	],
	entry_points={
		"console_scripts": [
			"gmc=gmc.__main__:main"
		]
	},
	package_data={
		"gmc.zzz": ["*.smk.py"]
	},
	include_package_data=True,
	data_files=[
		("gmc/etc", glob.glob("gmc/etc/*.*"))
	]
)

