from .mot import MOT
from .reading import smiles2graph, draw
import fire

def start():
    fire.Fire(MAPG)

class MAPG:
    
    def build_mot(self, smiles='CO', symmetry=False, mol_output='mol.svg', mot_output='mot.svg'):
        print(smiles)
        if mol_output is not None:
            svg = draw(smiles)
            with open(mol_output, 'w') as f:
                f.write(svg)
        mot = MOT(smiles, symmetry)
        mot.build()
        if mot_output is not None:
            with open(mot_output, 'w') as f:
                f.write(svg)