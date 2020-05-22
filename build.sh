#!/bin/bash

install_prefix=/ei/software/cb/minos
expected_version=$(grep -o "version *= *[\".0-9]\+" setup.py | sed "s/[\" ]//g" | cut -f 2 -d =)


if [[ -z "$expected_version" ]]; then 
 echo "Cannot parse version from setup.py. Please make sure to set a version (int)."
 exit
elif [[ "$#" -ge 1 ]]; then
  if [[ "$1" == "$expected_version" ]]; then
    version=$1
    if [[ "$#" -ge 2 ]]; then
      install_prefix=$2
    fi
  else
    echo "Supplied version ($1) does not match expected version ($expected_version)."
    exit
  fi
else
  version=$expected_version
fi

echo "Installing to $install_prefix/$version"



#rm dist/*whl
python setup.py bdist_wheel
pip install --prefix=${install_prefix}/${version}/x86_64 -U dist/minos-${version}-*.whl
#pip install --install-option="--prefix=/ei/software/testing/minos/${version}/x86_64" -U dist/minos-${version}-*.whl
