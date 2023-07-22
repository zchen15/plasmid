# Installation

## Requirements

Plasmid requires a number of Python packages, notably:

### Scientific libraries
- [NumPy](https://numpy.org/), for basic numerical routines
- [SciPy](https://www.scipy.org/), for sequence optimization
- [Pandas](https://pandas.pydata.org/), for manipulating data tables
- [biopython](https://biopython.org/), for reading bioinformatic files
- [scikit-learn](https://scikit-learn.org/), for sequence clustering

### Plotting libraries
- [matplotlib](https://matplotlib.org/), for generating static plots
- [bokeh](https://docs.bokeh.org/), for generating interactive plots
- [sty](https://pypi.org/project/sty/), for colorizing text

### Bioinformatic libraries
- [parasail](https://pypi.org/project/parasail/), for pairwise sequence alignment
- [pyspoa](https://pypi.org/project/pyspoa/), for multi-sequence alignment
- [mappy](https://pypi.org/project/mappy/), for python bindings to minimap2

### Bioinformatic software
Plasmid contains wrappers functions to a number of popular sequence alignment tools. The following are recommended to be installed on the system and accessible via commandline.
[minimap2](https://github.com/lh3/minimap2), for sequence alignment with long reads
[bwa](https://github.com/lh3/bwa), for sequence alignment with short reads
[mmseq2](https://github.com/soedinglab/MMseqs2), for database search

Plasmid is supported on Linux, macOS and Windows on Python 3.8 to 3.10.

## Installing via PyPI

You can also [install Plasmid from PyPI](https://pypi.org/project/plasmid/) using pip:

```
pip install plasmid
```

## Installing from github

Finally, you can also install the latest development version of plasmid from [GitHub](http://github.com/zchen15/Plasmid):

```
git clone git@github.com:zchen15/plasmid.git
cd plasmid
pip3 install .
```

This is useful if there is some feature that you want to try, but we did not release it yet as a stable version. Although you might find some unpolished details, these development installations should work without problems.
If you find any, please open an issue in the [issue tracker](https://github.com/zchen15/plasmid/issues).

```{warning}
It is recommended that you
**never ever use sudo** with distutils, pip, setuptools and friends in Linux
because you might seriously break your system
\[[1](http://wiki.python.org/moin/CheeseShopTutorial#Distutils_Installation)\]\[[2](http://stackoverflow.com/questions/4314376/how-can-i-install-a-python-egg-file/4314446#comment4690673_4314446)\]\[[3](http://workaround.org/easy-install-debian)\]\[[4](http://matplotlib.1069221.n5.nabble.com/Why-is-pip-not-mentioned-in-the-Installation-Documentation-tp39779p39812.html)\].
Use [virtual environments](https://docs.python.org/3/library/venv.html) instead.
```

## Problems and suggestions

If for any reason you get an unexpected error message or an incorrect result,
or you want to let the developers know about your use case,
please open a new issue in the [issue tracker](https://github.com/zchen15/plasmid/issues)
and we will try to answer promptly.