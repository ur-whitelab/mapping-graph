import networkx as nx
from .reading import smiles2graph
from .bond_equiv import bond_equiv_classes
import matplotlib.pyplot as plt
import io

class MOT:
    def __init__(self, smiles,symmetry=False):
        self._G, self._LG = smiles2graph(smiles)
        self._graph = nx.Graph()
        self._graph.add_node(frozenset(self._G.nodes()), {'atom_number': len(self._G.nodes())})
        #do symmetry processing if necessary
        self._bond_classes = bond_equiv_classes(self._G, self._LG)
        for bond in self._LG:
            index = -1
            degenerate = False
            for i,c in enumerate(self._bond_classes):
                if bond in c:
                    if len(c) > 1:
                        degenerate = True
                    index = i
                    break
            assert index != -1
            self._LG.node[bond]['equiv_class'] = self._bond_classes[index]
            self._LG.node[bond]['degenerate'] = degenerate

    def build(self):
        #build the MOT
        root = self._graph.nodes()[0]
        for bond in self._LG.nodes():
            self._remove_bond(self._G, self._LG, root, bond)
        #count number of unique number atoms, which corresponds to layer number
        self._layers = len(set([d['atom_number'] for n, d in self._graph.nodes(data=True)]))
        self._nodes = [[] for _ in range(self._layers)]
        for n,d in self._graph.nodes_iter(data=True):
            layer_i = d['atom_number'] - 1
            self._nodes[layer_i].append(set(n))

    def _remove_bond(self, G, LG, parent, bond, recurse=True, symmetry=False):
        #check for end
        if len(G) == 0:
            return

        #remove bond, check two atoms that make it up
        for i in range(2):
            #copy the grpahs
            _G, _LG = G.copy(), LG.copy()

            j = (i + 1) % 2
            #remove equivalent bonds from the graph and line graph
            if symmetry and self._LG.node[bond]['degenerate']:
                for equiv_bond in self._LG.node[bond]['equiv_class']:
                    #make sure we don't try to remove the bond twice
                    if equiv_bond != bond:
                        #actually remove the other atom and bond
                        if equiv_bond[j] in _G:
                            _G.remove_node(equiv_bond[j])
                        _LG.remove_node(equiv_bond)

            #actually remove the other atom and bond
            if bond[j] in _G:
                _G.remove_node(bond[j])
            _LG.remove_node(bond)
            #check if the graph is still connected
            if nx.is_connected(_G):
                new_node = frozenset(_G.nodes())
                self._graph.add_node(new_node, {'atom_number': len(_G.nodes())})
                self._graph.add_edge(parent, new_node)
                #now remove all bonds
                for new_bond in _LG.nodes():
                    if recurse:
                        if symmetry:
                            #need to only remove one bond per equivalence class
                            if new_bond != _LG.node[new_bond]['equiv_class'][0]:
                                continue
                        self._remove_bond(_G, _LG, new_node, new_bond)

    def __getitem__(self, index):
        #return as set so that it can be added
        return self._nodes[index[0]][index[1]]

    def validate_path(self, path):
        pass
    
    def draw(self):
        pos=nx.graphviz_layout(self._graph, prog="twopi")
        fig = plt.figure()
        nx.draw(self._graph, pos)
        with io.BytesIO() as output:
            fig.savefig(output, format='jpg')
            fig.clf()
            plt.close('all')
            return output.getvalue()