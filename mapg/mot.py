from .reading import smiles2graph, chem_line_graph
from .equiv import equiv_classes, equiv_function

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import io

import networkx as nx
import pygraphviz
from networkx.drawing.nx_agraph import graphviz_layout


class MOT:
    def __init__(self, smiles, symmetry=True, tree=True):
        self._G, self._LG = smiles2graph(smiles)
        self._graph = nx.DiGraph()
        self._graph.add_node(frozenset(self._G.nodes()), {'atom_number': len(self._G.nodes()),
                                                           'label': 'root'})
        self._symmetry = symmetry
        self._tree = tree
    def build(self):
        #build the MOT
        root = self._graph.nodes()[0]
        #create starting graph, which is quotient graph
        equiv_fxn = lambda a, b: False
        if self._symmetry:
            equiv_fxn = equiv_function(equiv_classes(self._G, 'atom_type', 'bond'))

        qgraph = nx.algorithms.minors.quotient_graph(self._G, equiv_fxn)
        self._graph.node[root]['atom_number'] = self._remove_bond(qgraph, self._G, root, symmetry=self._symmetry, tree=self._tree)

    def _remove_bond(self, graph, mol, parent, recurse=True, symmetry=False, layer=0, tree=False):
        #check for end

        if len(graph) == 1:
            return layer
        for nbunch in graph:
            atom_string = ','.join([mol.node[n]['atom_type'] for n in nbunch])
            id_string = ','.join([str(n) for n in nbunch])
            if nbunch not in self._graph:
                self._graph.add_node(nbunch, {'atom_number': layer,
                                            'label': atom_string + '\n' + id_string})
            else:
                self._graph.node[nbunch]['atom_number'] = min(self._graph.node[nbunch]['atom_number'], layer)
        #now remove bond
        for bond in graph.edges():
            _G = nx.contracted_edge(graph, bond, False)
            #we need to relabel it with all of the component nodes
            new_node = set(bond[0])
            new_node |= set(bond[1])
            nx.relabel_nodes(_G, {bond[0]: frozenset(new_node)}, copy=False)
            if recurse:
                last_layer = self._remove_bond(_G, mol, nbunch, recurse, symmetry, layer + 1, tree)
            if tree:
                break
        return last_layer

    def prune_parents(self):
        '''Reduce the number of parents'''
        possible_children = self._graph.nodes()
        possible_children.sort(key=lambda c: self._graph.node[c]['atom_number'], reverse=True)
        last_layer = self._graph.nodes(data=True)[0][1]['atom_number']
        layer_i = [None for _ in range(last_layer + 1)]
        for i,n in enumerate(possible_children):
            if last_layer != self._graph.node[n]['atom_number']:
                last_layer = self._graph.node[n]['atom_number']
                layer_i[last_layer] = i

        for n,d in self._graph.nodes_iter(data=True):
            #only keep edges that add to node set
            atom_number = self._graph.node[n]['atom_number']
            if atom_number == 0:
                continue
            #get the next lowest layer
            atom_number -= 1
            while layer_i[atom_number] is None:
                atom_number -= 1
            index = layer_i[atom_number]

            required = set(n)
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

    def validate_path(self, path):
        pass

    def draw(self, format='svg'):
        pos = graphviz_layout(self._graph, prog='dot', args='')
        fig = plt.figure(figsize=(12, 8))
        if len(self._graph) <= 25:
            nx.draw(self._graph, pos, labels={n: d['label'] for n,d in self._graph.nodes_iter(data=True)},
                node_size=5000)
        else:
            nx.draw(self._graph, pos)
        with io.BytesIO() as output:
            fig.savefig(output, format=format)
            fig.clf()
            plt.close('all')
            return output.getvalue()