# Module for molview viewer functionality
# 
# science libs
import os
import tempfile

import molviewspec as mvs

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
        self.work_dir = tempfile.TemporaryDirectory(prefix='molview_')

    def format_view(self, structure):
        '''
        Apply the default representations to an nglview widget.

        Args:
            view: nglview.widget.NGLWidget object.

        Returns:
            The formatted nglview widget.
        '''
        # colorize the structure
        (
            structure.component(selector="polymer")
                    .representation(type="cartoon")
                    .color(custom=dict(molstar_use_default_coloring=True))
        )
        (
            structure.component(selector="ligand")
                .representation(type="ball_and_stick")
                .color(custom=dict(molstar_use_default_coloring=True))
        )
        return structure
    
    def draw(self):
        '''
        Draw the record with molviewspec.
        
        Returns:
            tuple of (viewer, builder, structure).
        '''
        # generate the molviewspec
        builder = mvs.create_builder()
        builder.canvas(background_color=self.background)
        
        # format the structure
        structure = builder.download(url="mystructure.cif").parse(format="mmcif").model_structure()
        structure = self.format_view(structure)
        
        # write the structure to file
        filename = f'{self.work_dir.name}/record.cif'
        self.rec.save(filename)

        # render it with molstar
        with open(filename, "r") as f:
            cif_data = f.read()
        viewer = builder.molstar_notebook(data={"mystructure.cif": cif_data})
        
        self.builder = builder
        self.structure = structure
        self.viewer = viewer
        return viewer
        
    def draw_from_url(self, url:str, file_format:str):
        # generate the molviewspec
        builder = mvs.create_builder()
        builder.canvas(background_color=self.background)

        # fetch and format the structure
        structure = builder.download(url=url).parse(format=file_format).model_structure()
        structure = self.format_view(structure)
        
        viewer = builder.get_state()
        
        self.builder = builder
        self.structure = structure
        self.viewer = viewer
        return viewer
        