# Plasmid - A python library for gene editing and annotation
`Plasmid` is a library for reading, annotating, and editing genbank records in Jupyter Notebooks or ipython shell. This project was inspired by [ApE](https://jorgensen.biology.utah.edu/wayned/ape/) and [pandas](https://pandas.pydata.org/). The goal of the project was to address the lack of cross platform support provided by `ApE` and make the functionalities of `Biopython` more user friendly. This project is meant to serve as an open source interface for manipulating genbank records, querying genes, and manipulating NGS data in a programmable and user friendly manner that is cross compatible across all operating systems.

The functionalities of this library are organized into the following modules

`Plasmid` contains functions for reading, generating, and manipulating genbank records

`Graphic` contains functions for text colorization and plot generation

`Aligner` contains functions for sequence alignment and search

`Clust` contains functions for clustering sequences

`Designer` contains functions for generation of primers for extension PCR, Gibson Assembly, and Golden Gate Assembly.

## Installation
Run the following commands to clone and install this library. This package will be published to pypi once testing is complete.

```
git clone https://github.com/zchen15/plasmid
cd plasmid
pip3 install .
```

## Usage 
The `notebooks/` directory contains example jupyter notebooks that display some of functionality provided by `Plasmid`.

[genbank_editing.ipynb](https://zchen15.github.io/pages/misc/genbank_editing.html) shows how the `plasmid` module can be used to view genbank records, search for genes, and add new gene annotations.

[primer_design.ipynb](https://zchen15.github.io/pages/misc/primer_design.html) shows how the `designer` module can be used to generate primers for extension PCR, gibson assembly, and golden gate assembly.

[clustering.ipynb](https://zchen15.github.io/pages/misc/clustering.html) shows how the `aligner` and `cluster` modules can be used to align NGS reads, perform sequence clustering, multi-alignment, and database searches.

More detailed documentation will be provided later.

## Issues
If for any reason you get an unexpected error message or an incorrect result, or you want to let the developers know about your use case, please open a new issue in the [issue tracker](https://github.com/zchen15/plasmid/issues) and we will try to answer promptly.
