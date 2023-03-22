# Plasmid: A python library for gene editing and annotation
`Plasmid` is a python library for reading, annotating, and editing genbank records in Jupyter Notebooks or ipython shell. This project was inspired by [ApE](https://jorgensen.biology.utah.edu/wayned/ape/) and [pandas](https://pandas.pydata.org/). The goal of the project was address the lack of cross platform support provided by `ApE` and provide a more user friendly open source interface for manipulating genbank records in a programmable manner.

The functionalities of this library are organized into the following modules

`Plasmid` contains functions for reading, generating, and manipulating genbank records

`Graphic` contains functions for text colorization and plot generation

`Aligner` contains functions for sequence alignment and search

`Designer` contains functions for generation of primers for extension PCR, Gibson Assembly, and Golden Gate Assembly

## Installation
Run the following commands to clone and install this library. This package will be published to pypi once testing is complete.

```
git clone https://github.com/zchen15/plasmid
cd plasmid
pip3 install .
```

## Usage 
The `demos/` directory contains example jupyter notebooks that display some of functionality provided by `Plasmid`. More detailed documentation will be provided later.

## Issues
If you experience any issues with the code, please post them on the issues section along with the log file. I will monitor this periodically and try to fix issues as they arise.
