import networkx as nx
from .reading import smiles2graph

class MOT:
    def __init__(self, smiles):
        self._G, self._LG = smiles2graph(smiles)
        self._graph = nx.Graph()
        self._graph.add_node(frozenset(self._G.nodes()), {'atom_number': len(self._G.nodes())})

    def build(self):
        # do something to get MOT
        root = self._graph.nodes()[0]
        for bond in self._LG.nodes():
            self._remove_bond(self._G, self._LG, root, bond)
        #count number of unique number atoms, which corresponds to layer number
        self._layers = len(set([d['atom_number'] for n, d in self._graph.nodes(data=True)]))
        self._nodes = [[] for _ in range(self._layers)]
        for n,d in self._graph.nodes_iter(data=True):
            layer_i = d['atom_number'] - 1
            self._nodes[layer_i].append(set(n))
    def _remove_bond(self, G, LG, parent, bond, recurse=True):
        #check for end
        if len(G) == 0:
            return

        #remove bond, check two atoms that make it up
        for i in range(2):
            #copy the grpahs
            _G, _LG = G.copy(), LG.copy()
            #actually remove the other atom and bond
            _G.remove_node(bond[(i + 1) % 2])
            _LG.remove_node(bond)
            #check if the graph is still connected
            if nx.is_connected(_G):
                new_node = frozenset(_G.nodes())
                self._graph.add_node(new_node, {'atom_number': len(_G.nodes())})
                self._graph.add_edge(parent, new_node)
                #now remove all bonds
                for new_bond in _LG.nodes():
                    if recurse:
                        self._remove_bond(_G, _LG, new_node, new_bond)

    def __getitem__(self, index):
        #return as set so that it can be added
        return self._nodes[index[0]][index[1]]

    def validate_path(self, path):
        pass