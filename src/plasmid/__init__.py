#!/usr/bin/env python
# This is a python library for editing and plotting genbank files

# full in other libraries
from .genbank.plasmid import Plasmid
from .genbank.fileIO import fileIO
from .pcr.designer import Designer
from .pcr.oligos import Oligos
from .aligners.aligner import Aligner
from .cluster import Clust
from .graphic import Graphic

__version__ = "2026.11.3"
__author__ = "Zhewei Chen"
__email__ = "zchen15@github.com"
__name__ = "plasmid"
__description__ = 'Plasmid: A python package for viewing, editing, and assembling DNA'

