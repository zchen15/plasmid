#!/bin/bash

#mkdir sphinx
#cd sphinx
#sphinx-quickstart
#cd ..
#sphinx-apidoc -o sphinx .

# make output docs directory
mkdir docs
cd sphinx
make html
cp index.html ../docs/

