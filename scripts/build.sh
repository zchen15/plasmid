#!/bin/bash

# Build the sphinx documentation
echo installing latest version
pip3 install .

echo generating html docs with sphinx
mkdir docs

echo generating new docs
cd sphinx
make html
cp index.html ../docs/

cd ..

# Garbage collection
git checkout --orphan newbranch
git add -A
git commit -m 'cleanup'
git branch -D main
git branch -m main
git push -f origin main
git gc --aggressive --prune=all
