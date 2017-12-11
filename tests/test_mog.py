import networkx as nx
from mapg import *

def test_first_layer_is_atoms():
    mog = MOG('CO', symmetry=False)
    mol_graph = smiles2graph('CO')
    atoms = mog[0]
    for a in mol_graph:
        assert frozenset([a]) in mog[0]

    #now check that we get the number of equivalence atoms in layer
    mog = MOG('CO', symmetry=True)
    assert len(mog[0]) == 4

def test_layer_structure():
    mog = MOG('CCC', symmetry=True)
    for layer in mog:
        assert len(layer) > 0

    assert len(layer[-1]) == 1

def test_known_mog():
    mog = MOG('CO', symmetry=True)
    assert len(mog) == 4
    assert len(mog[0]) == 4
    assert len(mog[1]) == 3
    assert len(mog[2]) == 2
    assert len(mog[3]) == 1

