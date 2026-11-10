# Module for nglview viewer functionality

# science libs
import tempfile

import nglview
import mdtraj


def get_viewer_ngl(rec) -> 'nglview.widget.NGLWidget':
    '''
    Get an nglview object from the AtomArray

    Args:
        rec: Record object

    Returns:
        nglview.widget.NGLWidget object
    '''
    work_dir = tempfile.TemporaryDirectory(prefix='nglview_')
    filename = f'{work_dir.name}/record.pdb'
    rec.save(filename)
    view = nglview.show_file(filename)

    # default representations
    view.clear_representations()
    view.add_representation('cartoon', selection='protein', color='sstruc')
    view.add_representation('ball+stick', selection='ligand', color='element')
    view.camera = 'orthographic'
    view.background = 'gray'
    return view


def get_poses_ngl(rec) -> 'nglview.widget.NGLWidget':
    '''
    Get number of poses for nglview object from the AtomStack

    Args:
        rec: Record object

    Returns:
        nglview.widget.NGLWidget object
    '''
    work_dir = tempfile.TemporaryDirectory(prefix='nglview_')
    filename = f'{work_dir.name}/poses.pdb'
    rec.save(filename)

    traj = mdtraj.load(filename)
    view = nglview.show_mdtraj(traj)

    # default representations
    view.clear_representations()
    view.add_representation('cartoon', selection='protein', color='sstruc')
    view.add_representation('ball+stick', selection='ligand', color='element')
    view.camera = 'orthographic'
    view.background = 'gray'
    return view
