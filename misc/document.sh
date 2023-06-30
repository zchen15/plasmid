#!/bin/bash

mkdir docs
cd docs
sphinx-quickstart
cd ..
sphinx-apidoc -o docs .

