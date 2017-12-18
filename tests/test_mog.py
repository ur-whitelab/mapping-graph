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
    assert len(mog[-1]) == 1

def test_mog_isdag():
    mog = MOG('CO')
    assert nx.algorithms.dag.is_directed_acyclic_graph(mog.graph)

    mog = MOG('CC1=CC=CC=C1')
    assert nx.algorithms.dag.is_directed_acyclic_graph(mog.graph)

    mog = MOG('CCC')
    assert nx.algorithms.dag.is_directed_acyclic_graph(mog.graph)

    mog = MOG('CC')
    assert nx.algorithms.dag.is_directed_acyclic_graph(mog.graph)

def test_known_mog():
    mog = MOG('CO', symmetry=True)
    assert len(mog) == 4
    assert len(mog[0]) == 4
    assert len(mog[1]) == 3
    assert len(mog[2]) == 2
    assert len(mog[3]) == 1

def test_drawing_mog():
    mog = MOG('C1=CC=CC=C1')
    svg = mog.draw()

def test_path_matrix():
    mog = MOG('CO', symmetry=True)
    #validate paths are same size
    for p in mog.path_matrix:
        assert len(p) == len(mog.path_matrix[0])

def test_timeout():
    try:
        mog = MOG('CC1=CC=CC=C1', symmetry=False, build_timeout=1)
        assert False, 'Did not throw timeout'
    except MOG.TimeoutError:
        pass

