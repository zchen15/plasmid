# Module for py3Dmol viewer functionality

# science libs
import tempfile
import py3Dmol
import urllib

# custom libraries
from .sequences import get_sequence, get_ligands

class Viewer:
    '''
    This class holds functions pertaining to drawing structures with py3Dmol.

    Attributes:
        verbose: print run settings when True.
        vdw: show a van der Waals surface.
        surface_opacity: opacity of the van der Waals surface.
        surface_color: color of the van der Waals surface.
    '''
    verbose = False
    vdw = False
    height = 500
    width = 800
    surface_opacity = 0.7
    surface_color = 'white'
    background = 'black'

    def __init__(self, rec, **kwargs):
        '''
        Initialize the Viewer object and draw the record.
        
        Args:
            rec: Record object
            kwargs: settings that override the class defaults
        '''
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        self.rec = rec
        self.work_dir = tempfile.TemporaryDirectory(prefix='py3Dmol_')
        
    def format_view(self, view):
        '''
        Apply the default representations to an nglview widget.
        
        Args:
            view: nglview.widget.NGLWidget object.
        
        Returns:
            The formatted nglview widget.
        '''
        # colorize protein as ribbons
        df = get_sequence(self.rec).drop_duplicates()
        cols = ['chain_id', 'res_id']
        for chain_id, res_id in df[cols].values:
            protein = {'resi': f'{min(res_id.positions)}-{max(res_id.positions)}'}
            view.setStyle(protein, {'cartoon':{}})
            
            # show van der Waals surface
            if self.vdw:
                surface = {'opacity': self.surface_opacity,
                           'color': self.surface_color}
                view.addSurface(py3Dmol.VDW, surface, protein)
        
        # colorize ligands as sticks
        df = get_ligands(self.rec)
        for res_name in df['res_name'].values:
            ligand = {'resn': res_name}
            view.setStyle(ligand, {
                'stick': {'radius': 0.15, 'colorscheme': 'Jmol'}, 
                'sphere': {'scale': 0.25, 'colorscheme': 'Jmol'},
            })
        return view
        
    def draw(self):
        '''
        Draw the record with py3Dmol.

        Returns:
            py3Dmol.view object.
        '''
        # initialize the viewer object
        view = py3Dmol.view(width=self.width, height=self.height)
        view.setBackgroundColor(self.background)
        
        # load the data
        filename = f'{self.work_dir.name}/record.cif'
        self.rec.save(filename)
        with open(filename, 'r') as f:
            cif_data = f.read()
        view.addModel(cif_data, "mmcif")
        
        # format the view
        view = self.format_view(view)
        return view

    def draw_from_url(self, url:str, file_format:str):
        # initialize the viewer object
        view = py3Dmol.view(width=self.width, height=self.height)
        view.setBackgroundColor(self.background)

        # load the data
        with urllib.request.urlopen(url) as response:
            cif_data = response.read().decode('utf-8')
        view.addModel(cif_data, file_format)
        
        # format the view
        view = self.format_view(view)
        return view
