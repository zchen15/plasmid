#!/bin/bash

rm -rf dist build */*.egg-info *.egg-info
pip3 install .
pytest
