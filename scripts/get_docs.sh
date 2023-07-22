#!/bin/bash

#echo initializing sphinx docs
#mkdir sphinx
#cd sphinx
#sphinx-quickstart
#cd ..
#sphinx-apidoc -o sphinx src
# make output docs directory

echo generating html docs with sphinx
mkdir docs
cd sphinx
make html
cp index.html ../docs/
