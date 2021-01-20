#!/bin/bash
install_path=$1

export CPC_HOME=$install_path/CPC2-beta
export PATH=$PATH:$CPC_HOME/bin


cd $(dirname $CPC_HOME)
# http://cpc2.cbi.pku.edu.cn/download.php
wget http://cpc2.cbi.pku.edu.cn/data/CPC2-beta.tar.gz
tar xvzf CPC2-beta.tar.gz
rm CPC2-beta.tar.gz
cd CPC2-beta/libs/libsvm
tar xvzf libsvm-3.18.tar.gz
cd libsvm-3.18
make clean && make
cd $CPC_HOME/bin

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
