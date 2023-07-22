#!/bin/bash

#echo initializing sphinx docs
#mkdir sphinx
#cd sphinx
#sphinx-quickstart
#cd ..
#sphinx-apidoc -o sphinx src
# make output docs directory
echo installing latest version
pip3 install .

echo generating html docs with sphinx
mkdir docs

#echo clearing old docs
#rm -rf docs/html
#rm -rf docs/doctrees

echo generating new docs
cd sphinx
make html
cp index.html ../docs/
