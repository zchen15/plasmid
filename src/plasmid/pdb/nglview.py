# Module for nglview viewer functionality

import tempfile

from ..utils import library_safeload

mdtraj = library_safeload("mdtraj")
nglview = library_safeload("nglview")

from .molstar import feature_representation


class Viewer:
    '''
    This class holds functions pertaining to drawing structures with nglview.

    Attributes:
        verbose: print run settings when True.
        poses: draw a multi-model AtomStack as a trajectory.
        cartoon: show protein cartoons.
        ligand: show ligands as ball and stick.
        camera: camera type.
        background: background color.
    '''
    verbose = False
    camera = 'orthographic'
    background = 'black'

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
        self.work_dir = tempfile.TemporaryDirectory(prefix='nglview_')
        
    def format_view(self, view):
        '''
        Apply the default representations and features to an nglview widget.

        Args:
            view: nglview.widget.NGLWidget object.

        Returns:
            The formatted nglview widget.
        '''
        view.clear_representations()
        view.add_representation('cartoon', selection='protein')
        
        view.add_representation('basepair', selection='nucleic', color='red')
        view.add_representation('cartoon', selection='nucleic')
        #view.add_representation('basepair', selection='nucleic')
        
        view.add_representation('ball+stick', selection='ligand', color='element')

        # apply feature representations and colors
        if self.rec.features is not None and len(self.rec.features) > 0:
            if self.verbose:
                print(f"applying the following representations:\n{self.rec.features}")
            for _, feature in self.rec.features.iterrows():
                if len(feature["atom_id"]) == 0:
                    continue
                color = feature.get("color")
                if not isinstance(color, str) or not color.strip():
                    color = None
                kwargs = {}
                if color is not None:
                    kwargs["color"] = color.strip()
                view.add_representation(
                    feature_representation(feature),
                    selection=list(feature["atom_id"].positions),
                    **kwargs,
                )

        view.camera = self.camera
        view.background = self.background
        return view

    def draw(self):
        '''
        Draw a single structure from the record.

        Returns:
            nglview.widget.NGLWidget object.
        '''
        filename = f'{self.work_dir.name}/record.pdb'
        self.rec.save(filename)
        view = nglview.show_file(filename)
        view = self.format_view(view)
        return view

    def draw_from_url(self, url:str):
        '''
        Draw a single structure from the record.

        Returns:
            nglview.widget.NGLWidget object.
        '''
        view = nglview.show_url(url)
        view = self.format_view(view)
        return view
    
    def draw_poses(self):
        '''
        Draw a multi-model record as a trajectory.

        Returns:
            nglview.widget.NGLWidget object.
        '''
        filename = f'{self.work_dir.name}/poses.pdb'
        self.rec.save(filename)
        traj = mdtraj.load(filename)
        view = nglview.show_mdtraj(traj)
        view = self.format_view(view)
        return view
