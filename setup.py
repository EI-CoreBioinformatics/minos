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

name="minos"
version = "1.7.2"

if sys.version_info.major != 3:
	raise EnvironmentError("""minos is a python module that requires python3, and is not compatible with python2. Also, it is now 2020 and support for 2.x has ceased.""")

requirements = [line.rstrip() for line in open("requirements.txt", "rt")]

setup(
	name=name,
	python_requires=">=3.6",
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
        install_requires=requirements,
	entry_points={
		"console_scripts": [
			"minos=minos.__main__:main",
			"generate_metrics=minos.scripts.generate_metrics:main",
			"parse_mikado_gff=minos.scripts.parse_mikado_gff:main",
			"protein_completeness=minos.scripts.protein_completeness:main",
			"collapse_metrics=minos.scripts.collapse_metrics:main",
			"validate_gff3=minos.scripts.validate_gff3:main",
			"create_release_gff3=minos.scripts.create_release_gff:main",
			"sanity_check=minos.scripts.sanity_check:main",
			"parse_mikado_stats=minos.scripts.parse_mikado_stats:main",
			"analyse_busco=minos.scripts.analyse_busco:main",
			"install_cpc2=minos.scripts.install_cpc2:main"
		]
	},
	package_data={
		"minos.zzz": ["*.smk"],
		"minos.etc": ["*.json", "*.yaml", "*.def", "test/*"],
		"minos.dependencies": ["*.sh", "container_defs/*"]
	},
	include_package_data=True,
)
