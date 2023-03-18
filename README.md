
## Introduction
The following is a python library for editing plasmids in jupyter notebook or ipython shell. This library was inspired by [ApE](https://jorgensen.biology.utah.edu/wayned/ape/) and [pandas](https://pandas.pydata.org/). Genbank information is stored in an object class called plasmid. The plasmid object can be manipulated in similar manner to a pandas dataframe. Useful functions for gene editing have been supplemented in this class.

## Dependencies
This library uses [biopython](https://biopython.org/), [dna_features_viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer), and [pyspoa](https://github.com/nanoporetech/pyspoa).
```
pip install biopython             # for reading genbank files
pip install pyspoa                # for multi-alignment
```

## Installation
Clone this repository and link `gene_editor.py` to the current directory where you will run your jupyter notebook.

## Quick start
```
import gene_editor as pyge        # loads the library

plasmid = pyge.read('seq.gb')     # loads genbank data into the plasmid dataframe
print(plasmid)                    # prints a table containing the layout of the plasmid
df = plasmid.get_dataframe()      # obtain a pandas dataframe of the genetic features on the plasmid

graphic_record = plasmid.linear_record() # obtain a linear graphic record for plotting with dna_features_viewer
plasmid.linear_record().plot()           # obtain a linear graphic record and plot a linear map of plasmid
plasmid.circular_record().plot()         # obtain a circular graphic record and plot a circular map of the plasmid

plasmid.to_genbank('output.gb')  # writes the plasmid dataframe to a genbank file 
```

More examples can be found in the demos folder.

## Issues
Please submit issues on the issues tab of this repository.

## Citation
If you use this library for scientific publication, please also cite [dna_features_viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer)


