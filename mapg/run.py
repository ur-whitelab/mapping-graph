from .mot import MOT
from .reading import smiles2graph, draw
from .bond_equiv import bond_equiv_classes
import fire

def start():
    fire.Fire({
        'build-mot': build_mot
    })

def build_mot(smiles='CO', symmetry=True, mol_output='mol.svg', mot_output='mot.svg'):
    if mol_output is not None:
        bec = bond_equiv_classes(*smiles2graph(smiles))
        svg = draw(smiles, bec, True)
        with open(mol_output, 'w') as f:
            f.write(svg)
    mot = MOT(smiles, symmetry)
    mot.build()
    mot.prune_parents()
    mot.prune_nodes()
    if mot_output is not None:
        svg = mot.draw()
        with open(mot_output, 'wb') as f:
            f.write(svg)