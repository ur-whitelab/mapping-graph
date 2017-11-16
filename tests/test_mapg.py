import mapg

def test_smiles2graph():
    '''Tests a basic conversion from smiles into a graph'''
    graph = mapg.smiles2graph('CO')
    assert(True)