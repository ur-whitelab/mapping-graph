from .reading import smiles2graph, chem_line_graph
from .equiv import equiv_classes

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
        self._bond_classes = equiv_classes(self._G, self._LG)
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
                            _G.remove_node(equiv_bond[i])
                            #pass
                        _LG.remove_node(equiv_bond)

            #actually remove the other atom and bond
            if bond[j] in _G:
                _G.remove_node(bond[j])
            _LG.remove_node(bond)
            #check if the graph is still connected
            for _SG in nx.connected_component_subgraphs(_G):
                if not (bond[0] in _SG or bond[1] in _SG):
                    continue
                _SLG = chem_line_graph(_SG)
                #copy over attributes
                for n in _SLG.nodes():
                    _SLG.node[n]['equiv_class'] = _LG.node[n]['equiv_class']
                    _SLG.node[n]['degenerate'] = _LG.node[n]['degenerate']
                new_node = frozenset(_SG.nodes())
                #new node will have the list of atoms that make it up
                #and the number of atoms

                #FYI node may already be there
                atom_string = ','.join([d['atom_type'] for _, d in _SG.nodes_iter(data=True)])
                id_string = ','.join([str(n) for n in _SG])
                self._graph.add_node(new_node, {'atom_number': len(_SG.nodes()),
                                                'label': atom_string + '\n' + id_string,
                                                'layer': layer})
                #add parent
                do_add = False
                for n in new_node:
                    do_add |= n in parent
                #if do_add:
                #self._graph.add_edge(parent, new_node)
                #now remove all bonds
                for new_bond in _SLG.nodes():
                    if recurse:
                        if symmetry:
                            #need to only remove one bond per equivalence class
                            if _SLG.node[new_bond]['degenerate']:
                                continue
                        self._remove_bond(_SG, _LG, new_node, new_bond, symmetry=symmetry, layer = layer + 1)

    def prune_parents(self):
        '''Reduce the number of parents'''
        possible_children = self._graph.nodes()
        possible_children.sort(key=lambda c: self._graph.node[c]['atom_number'], reverse=True)
        layer_i = [None for _ in range(len(self._G) + 1)]
        last_layer = len(self._G)
        for i,n in enumerate(possible_children):
            print(i, last_layer, self._graph.node[n]['atom_number'])
            if last_layer != self._graph.node[n]['atom_number']:
                last_layer = self._graph.node[n]['atom_number']
                layer_i[last_layer] = i
        print(layer_i)

        for n,d in self._graph.nodes_iter(data=True):
            #only keep edges that add to node set
            atom_number = self._graph.node[n]['atom_number']
            if atom_number == 1:
                continue
            #get the next lowest layer
            atom_number -= 1
            while layer_i[atom_number] is None:
                atom_number -= 1
            index = layer_i[atom_number]

            required = set(n)
            print(index, self._graph.node[n]['atom_number'])
            print(n, possible_children[index:])
            for c in possible_children[index:]:
                if not required.isdisjoint(c):
                    required -= c
                    self._graph.add_edge(n, c)
            #assert len(required ) == 0


    def prune_nodes(self):
        '''Remove nodes that have one parent and one child'''
        success = True
        while success:
            to_del = []
            for n in self._graph:
                if len(self._graph.in_edges(n)) == 1 and len(self._graph.out_edges(n)) == 1:
                    #add the edge now
                    self._graph.add_edge(list(self._graph.in_edges(n))[0][0],
                                        list(self._graph.out_edges(n))[0][1])
                    to_del.append(n)
            self._graph.remove_nodes_from(to_del)
            success = len(to_del) > 0

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