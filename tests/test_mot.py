import networkx as nx
from mapg import *

def test_remove_bond_mot():
    mot = MOT('CO')
    #try removing one bond from methanol
    #should have 3x CH2-OH, 1x CH3-O
    for bond in mot._LG.nodes():
        mot._remove_bond(mot._G, mot._LG, mot._graph.nodes()[0], bond, False)
    assert len(mot._graph) == 3 + 1 + 1

def test_remove_bond_mot_symm():
    mot = MOT('CO', symmetry=True)
    #try removing one bond from methanol
    #should have 1x CH2-OH, 1x CH3-O, and root
    for bond in mot._LG.nodes():
        mot._remove_bond(mot._G, mot._LG, mot._graph.nodes()[0], bond, False, True)
    assert len(mot._graph) == 1 + 1 + 1

def test_exhaustive_mot():
    '''
    Test that all CG mappings are accounted for in the MOT
    '''

def test_path_validator():
    '''
    Test that the validator correctly disallows certain node selections
    '''