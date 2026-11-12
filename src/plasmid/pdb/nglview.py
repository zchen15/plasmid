# Module for nglview viewer functionality

# science libs
import tempfile

import nglview
import mdtraj


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
        traj: cached mdtraj trajectory loaded for pose drawing.
    '''
    verbose = False
    poses = False
    cartoon = True
    ligand = True
    camera = 'orthographic'
    background = 'gray'
    traj = None

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
        self.view = self.draw()
        # keep the temporary directory alive for as long as the widget exists
        self.view._work_dir = self.work_dir

    def draw(self):
        '''
        Draw the record with nglview.

        Returns:
            nglview.widget.NGLWidget object.
        '''
        if self.poses:
            return self.draw_poses()
        return self.draw_structure()

    def format_view(self, view):
        '''
        Apply the default representations to an nglview widget.

        Args:
            view: nglview.widget.NGLWidget object.

        Returns:
            The formatted nglview widget.
        '''
        view.clear_representations()
        if self.cartoon:
            view.add_representation('cartoon', selection='protein', color='sstruc')
        if self.ligand:
            view.add_representation('ball+stick', selection='ligand', color='element')
        view.camera = self.camera
        view.background = self.background
        return view

    def draw_structure(self):
        '''
        Draw a single structure from the record.

        Returns:
            nglview.widget.NGLWidget object.
        '''
        filename = f'{self.work_dir.name}/record.pdb'
        self.rec.save(filename)
        return self.format_view(nglview.show_file(filename))

    def draw_poses(self):
        '''
        Draw a multi-model record as a trajectory.

        Returns:
            nglview.widget.NGLWidget object.
        '''
        if self.traj is None:
            filename = f'{self.work_dir.name}/poses.pdb'
            self.rec.save(filename)
            self.traj = mdtraj.load(filename)
        return self.format_view(nglview.show_mdtraj(self.traj))
