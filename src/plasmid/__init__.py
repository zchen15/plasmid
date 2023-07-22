#!/usr/bin/env python
# This is a python library for editing and plotting genbank files

# full in other libraries
from .plasmid import Plasmid
from .plasmid import read_genbank
from .designer import Designer
from .oligos import Oligos
from .aligner import Aligner
from .cluster import Clust
from .fileIO import fileIO
from .graphic import Graphic

__version__ = "2023.7.21"
__author__ = "Zhewei Chen"
__email__ = "zchen15@github.com"
__name__ = "plasmid"
__description__ = 'Plasmid: A python package for viewing, editing, and assembling DNA'

