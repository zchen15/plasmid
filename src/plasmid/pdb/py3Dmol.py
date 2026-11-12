# Module for py3Dmol viewer functionality

# science libs
import tempfile

import py3Dmol


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
    vdw = True
    surface_opacity = 0.7
    surface_color = 'white'

    def __init__(self, rec, **kwargs):
        '''
        Initialize the Viewer object and draw the record.

        Args:
            rec: Record object.
            kwargs: settings that override the class defaults.

        '''
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        self.rec = rec
        self.work_dir = tempfile.TemporaryDirectory(prefix='py3Dmol_')
        self.view = self.draw()

    def draw(self):
        '''
        Draw the record with py3Dmol.

        Returns:
            py3Dmol.view object.
        '''
        filename = f'{self.work_dir.name}/record.cif'
        self.rec.save(filename)

        with open(filename, 'r') as f:
            data = f.read()
        view = py3Dmol.view(data=data)

        # colorize protein as ribbons
        df = self.rec.get_sequence().drop_duplicates()
        cols = ['chain_id', 'res_id']
        for chain_id, res_id in df[cols].values:
            protein = {'resi': f'{min(res_id.positions)}-{max(res_id.positions)}'}
            view.setStyle(protein, {'cartoon': {'color': 'spectrum'}})

            # show van der Waals surface
            if self.vdw:
                surface = {'opacity': self.surface_opacity,
                           'color': self.surface_color}
                view.addSurface(py3Dmol.VDW, surface, protein)

        # colorize ligands as sticks
        df = self.rec.get_ligands()
        for res_name in df['res_name'].values:
            ligand = {'resn': res_name}
            view.setStyle(ligand, {'stick': {}})
        return view
