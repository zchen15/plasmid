# Module for molviewspec viewer functionality

# science libs
import tempfile

import molviewspec


def get_viewer_mvs(rec, height: int = 500, width: int = 700):
    '''
    Get molviewspec object from the AtomArray

    Args:
        rec: Record object
        height: height of the viewer
        width: width of the viewer

    Returns:
        molviewspec viewer object
    '''
    work_dir = tempfile.TemporaryDirectory(prefix='molviewspec_')
    filename = f'{work_dir.name}/record.cif'
    rec.save(filename)

    # read the data
    with open(filename, 'r') as f:
        cif_data = f.read()

    # generate the view
    builder = molviewspec.create_builder()
    structure = builder.download(url='mystructure.cif').parse(format='mmcif').model_structure()

    structure.component(selector='polymer')\
        .representation(type='cartoon')\
        .color(custom=dict(molstar_use_default_coloring=True))

    structure.component(selector='ligand')\
        .representation(type='ball_and_stick')\
            .color(custom=dict(molstar_use_default_coloring=True))

    # output the results
    viewer = builder.molstar_notebook(data={'mystructure.cif': cif_data}, width=width, height=height)
    return viewer, builder, structure
