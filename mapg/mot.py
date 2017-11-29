from .reading import smiles2graph
from .bond_equiv import bond_equiv_classes

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import io

import networkx as nx
import pygraphviz
from networkx.drawing.nx_agraph import graphviz_layout


class MOT:
    def __init__(self, smiles,symmetry=False):
        self._G, self._LG = smiles2graph(smiles)
        self._graph = nx.DiGraph()
        self._graph.add_node(frozenset(self._G.nodes()), {'atom_number': len(self._G.nodes()),
                                                           'label': 'root',
                                                           'layer': 0})
        #do symmetry processing if necessary
        self._bond_classes = bond_equiv_classes(self._G, self._LG)
        first_in_class = [True for _ in self._bond_classes]
        for bond in self._LG:
            index = -1
            degenerate = True
            for i,c in enumerate(self._bond_classes):
                if bond in c:
                    if first_in_class[i]:
                        degenerate = False
                        first_in_class[i] = False
                    index = i
                    break
            assert index != -1
            self._LG.node[bond]['equiv_class'] = self._bond_classes[index]
            self._LG.node[bond]['degenerate'] = degenerate #the chosen representative or not of bond class
        self._symmetry = symmetry
    def build(self):
        #build the MOT
        root = self._graph.nodes()[0]
        for bond in self._LG.nodes():
            if self._symmetry:
                #need to only remove one bond per equivalence class
                if self._LG.node[bond]['degenerate']:
                    continue
            self._remove_bond(self._G, self._LG, root, bond, symmetry=self._symmetry)
        #count number of unique number atoms, which corresponds to layer number
        self._layers = len(set([d['atom_number'] for n, d in self._graph.nodes(data=True)]))
        self._nodes = [[] for _ in range(self._layers)]
        for n,d in self._graph.nodes_iter(data=True):
            self._nodes[d['layer']].append(set(n))

    def _remove_bond(self, G, LG, parent, bond, recurse=True, symmetry=False, layer = 0):
        #check for end
        if len(G) == 0:
            return

        print(bond, layer)
        #remove bond and form two subsets from each side
        #of the two atoms that make it up
        for i in range(2):
            #copy the grpahs
            _G, _LG = G.copy(), LG.copy()

            j = (i + 1) % 2
            #remove equivalent bonds from the graph and line graph
            if symmetry:
                for equiv_bond in self._LG.node[bond]['equiv_class']:
                    #make sure we don't try to remove the main bond
                    if equiv_bond != bond:
                        #remove all atoms and bonds here,
                        #since these are the degenerate bonds
                        if equiv_bond[j] in _G:
                            _G.remove_node(equiv_bond[j])
                        #make sure we don't remove atom i if it's part of
                        #the designated non-degenerate bond
                        if equiv_bond[i] in _G and not equiv_bond[i] == bond[i]:
                            #_G.remove_node(equiv_bond[i])
                            pass
                        _LG.remove_node(equiv_bond)

            #actually remove the other atom and bond
            if bond[j] in _G:
                _G.remove_node(bond[j])
            _LG.remove_node(bond)
            #check if the graph is still connected
            print(layer, self._LG.node[bond]['bond'], i, j, _G.nodes(data=True))
            if len(_G) > 0:
                _G = list(nx.connected_component_subgraphs(_G))[0]
            if len(_G) > 0 and nx.is_connected(_G):
                new_node = frozenset(_G.nodes())
                #new node will have the list of atoms that make it up
                #and the number of atoms

                #FYI node may already be there
                atom_string = ','.join([d['atom_type'] for _, d in _G.nodes_iter(data=True)])
                id_string = ','.join([str(n) for n in _G])
                self._graph.add_node(new_node, {'atom_number': len(_G.nodes()),
                                                'label': atom_string + '\n' + id_string,
                                                'layer': layer})
                #only add if we are the first to add it
                self._graph.add_edge(parent, new_node)
                #now remove all bonds
                for new_bond in _LG.nodes():
                    if recurse:
                        if symmetry:
                            #need to only remove one bond per equivalence class
                            if _LG.node[new_bond]['degenerate']:
                                continue
                        self._remove_bond(_G, _LG, new_node, new_bond, symmetry=symmetry, layer = layer + 1)

    def __getitem__(self, index):
        #return as set so that it can be added
        return self._nodes[index[0]][index[1]]

    def validate_path(self, path):
        pass

    def draw(self):
        pos = graphviz_layout(self._graph, prog='dot', args='')
        fig = plt.figure()
        nx.draw(self._graph, pos, labels={n: d['label'] for n,d in self._graph.nodes_iter(data=True)})
        with io.BytesIO() as output:
            fig.savefig(output, format='svg')
            fig.clf()
            plt.close('all')
            return output.getvalue()