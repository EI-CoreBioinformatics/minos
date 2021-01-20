#!/usr/bin/env python
import stat
import re
import sys
import os
import subprocess
import argparse
import glob

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



def main():

	ap = argparse.ArgumentParser()
	ap.add_argument("--hash", default="9a2b596", help="Hash or version tag to use for downloading CPC2 from github")
	ap.add_argument("--repo-name", default="CPC2", help="Name of the repository in GitHub to download from.")
	ap.add_argument("-u", "--username", default="biocoder", help="Username associated with the CPC2 repo on GitHub")
	ap.add_argument("install_path", default=".", type=str, nargs="?")
	args = ap.parse_args()

	args.install_path = os.path.abspath(args.install_path)

	args.hash = args.hash[:7]
	with cd(args.install_path):
		os.makedirs("CPC2", exist_ok=True)
		cmd = "wget -O CPC2.tar.gz --no-check-certificate --content-disposition https://github.com/{args.username}/{args.repo_name}/tarball/{args.hash}".format(
			args=args)
		print("Downloading CPC2...", end="", flush=True)
		subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
		print(" done.")
		cmd = "tar --strip-components=1 -C CPC2 -xvzf CPC2.tar.gz".format(args=args)
		print("Unpacking...", end="", flush=True)
		subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
		print(" done.")

		cmd = "rm -f CPC2.tar.gz".format(args=args)
		print("Removing tarball...", end="", flush=True)
		subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
		print(" done.")

		cmd = "cd CPC2/libs/libsvm/libsvm-3.18 && make clean && make".format(args=args)
		print("Unpacking and building libsvm...", end="", flush=True)
		subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
		print(" done.")

	print("Patching CPC2...", end="", flush=True)
	cpc2_dir = os.path.join(args.install_path, "CPC2")
	bin_dir = os.path.join(cpc2_dir, "bin")
	compress_py = os.path.join(bin_dir, "compress.py")
	code = open(compress_py).read()
	code = re.sub("\\bfile[(]", "open(", code)
	with open(compress_py, "w") as _out:
		print(code, end="", file=_out, flush=True)
	
	cpc2_py = os.path.join(bin_dir, "CPC2.py")
	code = open(cpc2_py).read()
	code = re.sub("commands", "subprocess", code)
	code = re.sub("\\bfile[(]", "open(", code)
	code = re.sub("triplet_got\.next\(\)", "next(triplet_got)", code)
	with open(cpc2_py, "w") as _out:
		for i, line in enumerate(code.split("\n")):
			print(line, file=_out, flush=True)
			if i == 355:
				print("\t\tos.system('mv ' + outfile  + '.txt ' + outfile)", file=_out, flush=True)
	
	seqio_py = os.path.join(bin_dir, "seqio.py")
	code = open(seqio_py).read()
	# code = re.sub("print \([^ $]+\)", "print(\1)", code)
	code = re.sub("print (?P<prt>[a-zA-Z0-9_()]+)", "print(\g<prt>)", code)
	with open(seqio_py, "w") as _out:
		print(code, end="", file=_out, flush=True)		

	print(" done.")
	print("Changing file permissions...", end="", flush=True)
	for f in glob.glob(os.path.join(bin_dir, "*.py")):
		os.chmod(f, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
	print(" done.")

	print("CPC_HOME={}".format(cpc2_dir))
	print("PATH={}:{}/bin".format(os.environ.get("PATH"), cpc2_dir))

	pass


if __name__ == "__main__":
	main()

"""
export CPC_HOME=$install_path/CPC2-beta
export PATH=$PATH:$CPC_HOME/bin

# apply patches for Python3 compatibility
# patch compress.py
sed -i "s/\bfile(/open(/g" compress.py
# patch CPC2.py
sed -i "s/commands/subprocess/g" CPC2.py
sed -i "s/\bfile(/open(/g" CPC2.py
sed -i "s/triplet_got.next()/next(triplet_got)/g" CPC2.py

head -n 355 CPC2.py >> CPC2.py.1
echo $'\t'$'\t'"os.system('mv ' + outfile + '.txt ' + outfile)" >> CPC2.py.1
tail -n +356 CPC2.py >> CPC2.py.1
mv CPC2.py.1 CPC2.py
# patch seqio.py
sed -i "s/print \([^ $]\+\)/print(\1)/g" seqio.py
chmod 775 *.py

python CPC2.py -h
"""
