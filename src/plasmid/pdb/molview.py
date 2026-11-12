# Module for molview viewer functionality

# science libs
import os
import tempfile

import molviewspec


class Viewer:
    '''
    This class holds functions pertaining to drawing structures with molviewspec.

    Attributes:
        verbose: print run settings when True.
        height: height of the viewer.
        width: width of the viewer.
    '''
    verbose = False
    height = 500
    width = 700

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
        self.work_dir = tempfile.TemporaryDirectory(prefix='molview_')
        self.viewer, self.builder, self.structure = self.draw()

    def draw(self):
        '''
        Draw the record with molviewspec.

        Returns:
            tuple of (viewer, builder, structure).
        '''
        filename = f'{self.work_dir.name}/crystal.cif'
        self.rec.save(filename)

        # read the data
        with open(filename, 'r') as f:
            cif_data = f.read()

        # generate the view
        builder = molviewspec.create_builder()
        structure = builder.download(url='crystal.cif').parse(format='mmcif').model_structure()

        structure.component(selector='polymer')\
            .representation(type='cartoon')\
            .color(custom=dict(molstar_use_default_coloring=True))

        structure.component(selector='ligand')\
            .representation(type='ball_and_stick')\
            .color(custom=dict(molstar_use_default_coloring=True))

        # output the results
        viewer = builder.molstar_notebook(data={name: cif_data},
                                          width=self.width, height=self.height)
        return viewer, builder, structure
