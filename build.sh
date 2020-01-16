#!/bin/bash

version=$1

#rm dist/*whl
python setup.py bdist_wheel
pip install --prefix=/ei/software/testing/gmc/${version}/x86_64 -U dist/gmc-${version}-*.whl
