import mapg

def test_smiles2graph():
    '''Tests a basic conversion from smiles into a graph'''
    graph, line_graph = mapg.smiles2graph('CO')

    #check that the atom number is correct
    assert(len(graph.nodes()) == 6)

    # get index of carbon atom
    c_index, o_index = -1, -1
    for n,d in graph.nodes(data=True):
        if d['atomType'] == 'C':
            c_index = n
        elif d['atomType'] == 'O':
            o_index = n
    assert(graph.has_edge(c_index, o_index))

def test_draw():
    '''Smoke test'''
    svg = mapg.draw('CCCO')
    svg2 = mapg.draw('CCCO', [(1, 2), (3, 4)], True)

