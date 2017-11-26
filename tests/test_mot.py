import networkx as nx
from mapg import *

def test_remove_bond_mot():
    mot = MOT('CO')
    #try removing one bond from methanol
    #should have 3x CH2-OH, 1x CH3-O
    for bond in mot._LG.nodes():
        mot._remove_bond(mot._G, mot._LG, mot._graph.nodes()[0], bond, False)
    assert len(mot._graph) == 3 + 1 + 1

def test_partition_mot():
    '''
    Test that each layer is a valid partition of the atoms
    '''
    mot = MOT('CO')
    mot.build()
    layer = set()
    #addition is union, selection gets a layer
    for n in mot[1, :]:
        layer += n
    molgraph, _ = smiles2graph('CO')
    assert len(layer) == len(molgraph)



def test_exhaustive_mot():
    '''
    Test that all CG mappings are accounted for in the MOT
    '''

def test_path_validator():
    '''
    Test that the validator correctly disallows certain node selections
    '''