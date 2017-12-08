import mapg
import networkx as nx
import operator
def test_smiles2graph():
    '''Tests a basic conversion from smiles into a graph'''
    graph, line_graph = mapg.smiles2graph('CO')

    #check that the atom number is correct
    assert(len(graph.nodes()) == 6)

    # get index of carbon atom
    c_index, o_index = -1, -1
    for n,d in graph.nodes(data=True):
        if d['atom_type'] == 'C':
            c_index = n
        elif d['atom_type'] == 'O':
            o_index = n
    assert(graph.has_edge(c_index, o_index))

def test_draw():
    '''Smoke test'''
    svg = mapg.draw('CCCO')
    svg2 = mapg.draw('CCCO', [(1, 2), (3, 4)], True)


def test_bond_equiv_classes():
    '''Tests if the function correctly identifies
    equivalent bonds'''
    #check the number of equivalent bonds with no local or global symmetry
    graph1,line_graph1=mapg.smiles2graph('C(F)(Cl)Br')
    equiv_bond_groups1=mapg.equiv_classes(graph1,line_graph1)
    #Since there is no symmetry, all the four bonds will have their individual equivalence groups
    assert(len(equiv_bond_groups1)==4)

    #check number of bond classes in chloroacetone containing local symmetry
    graph2, line_graph2=mapg.smiles2graph('CC(=O)CCl')
    equiv_bond_groups2=mapg.equiv_classes(graph2,line_graph2)
    assert(len(equiv_bond_groups2)==6)

    #check number of bond classes in neopentane which displays global symmetry
    graph3, line_graph3=mapg.smiles2graph('CC(C)(C)C')
    equiv_bond_groups3=mapg.equiv_classes(graph3,line_graph3)
    assert(len(equiv_bond_groups3)==2)



