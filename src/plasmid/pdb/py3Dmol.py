# Module for py3Dmol viewer functionality

# science libs
import tempfile

import py3Dmol


def get_viewer_py3Dmol(rec, vdw: bool = True):
    '''
    Get a py3Dmol viewer object from the AtomArray

    Args:
        rec: Record object
        vdw: show vanderwaals surface

    Returns:
        py3Dmol.view object
    '''
    work_dir = tempfile.TemporaryDirectory(prefix='py3Dmol_')
    filename = f'{work_dir.name}/record.cif'
    rec.save(filename)

    with open(filename, 'r') as f:
        data = f.read()
    view = py3Dmol.view(data=data)

    # colorize protein as ribbons
    df = rec.get_sequence().drop_duplicates()
    cols = ['chain_id', 'res_id']
    for chain_id, res_id in df[cols].values:
        protein = {'resi': f'{min(res_id.positions)}-{max(res_id.positions)}'}
        view.setStyle(protein, {'cartoon': {'color': 'spectrum'}})

        # show vanderwaals
        if vdw:
            view.addSurface(py3Dmol.VDW, {'opacity': 0.7, 'color': 'white'}, protein)

    # colorize ligands as sticks
    df = rec.get_ligands()
    for res_name in df['res_name'].values:
        ligand = {'resn': res_name}
        view.setStyle(ligand, {'stick': {}})
    return view
