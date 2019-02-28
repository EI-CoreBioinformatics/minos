#!/bin/bash

rm dist/*whl
python setup.py bdist_wheel
pip install --prefix=/ei/software/testing/gmc/dev/x86_64 -U dist/gmc*.whl
