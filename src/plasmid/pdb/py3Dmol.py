# Module for py3Dmol viewer functionality

# science libs
import tempfile
import urllib.request

import molviewspec as mvs
import py3Dmol

from .molstar import feature_representation

# Mapping of molstar representation types to py3Dmol.js style names.
REPRESENTATION_STYLES = {
    "cartoon": "cartoon",
    "backbone": "cartoon",
    "ball_and_stick": "stick",
    "spacefill": "sphere",
    "line": "lines",
    "surface": "surface",
}


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
        Apply the default representations and features to a py3Dmol view.

        Args:
            view: py3Dmol.view object.

        Returns:
            The formatted py3Dmol view.
        '''
        # colorize protein as ribbons
        for method in (self.rec.get_amino_acids, self.rec.get_nucleic_acids):
            df = method()
            seen = set()
            for _, row in df.iterrows():
                chain_id = str(row['chain_id'])
                res_id = row['res_id']
                key = (chain_id, tuple(int(x) for x in res_id.positions))
                if key in seen:
                    continue
                seen.add(key)
                protein = {'resi': f'{min(res_id.positions)}-{max(res_id.positions)}'}
                view.setStyle(protein, {'cartoon':{}})
                
                # show van der Waals surface
                if self.vdw:
                    surface = {'opacity': self.surface_opacity,
                               'color': self.surface_color}
                    view.addSurface(py3Dmol.VDW, surface, protein)
        
        # colorize ligands as sticks
        df = self.rec.get_ligands()
        for res_name in df['res_name'].values:
            ligand = {'resn': res_name}
            view.setStyle(ligand, {
                'stick': {'radius': 0.15, 'colorscheme': 'Jmol'}, 
                'sphere': {'scale': 0.25, 'colorscheme': 'Jmol'},
            })
        
        # apply feature representations and colors
        if self.rec.features is not None and len(self.rec.features) > 0:
            if self.verbose:
                print(f"applying the following representations:\n{self.rec.features}")
            for _, feature in self.rec.features.iterrows():
                selectors = [
                    mvs.ComponentExpression(
                        auth_asym_id=str(feature["chain_id"]),
                        beg_auth_seq_id=int(s.start),
                        end_auth_seq_id=int(s.stop) - 1,
                    )
                    for s in feature["res_id"].get_slices()
                ]
                if len(selectors) == 0:
                    continue
                representation_type = feature_representation(feature)
                color = feature.get("color")
                if not isinstance(color, str) or not color.strip():
                    color = None

                for expr in selectors:
                    selection = {}
                    if expr.auth_asym_id is not None:
                        selection["chain"] = str(expr.auth_asym_id)
                    if expr.beg_auth_seq_id is not None and expr.end_auth_seq_id is not None:
                        selection["resi"] = f"{expr.beg_auth_seq_id}-{expr.end_auth_seq_id}"

                    if representation_type == "cartoon":
                        style = {"cartoon": {"color": color}} if color is not None else {"cartoon": {}}
                    elif representation_type == "ball_and_stick":
                        style = {
                            'stick': {'radius': 0.15, 'colorscheme': color or 'Jmol'},
                            'sphere': {'scale': 0.25, 'colorscheme': color or 'Jmol'},
                        }
                    else:
                        style = {REPRESENTATION_STYLES.get(representation_type, "cartoon"): {}}
                    view.setStyle(selection, style)
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
