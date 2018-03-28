from .reading import smiles2graph, chem_line_graph
from .equiv import equiv_classes, equiv_function
from .util import integer_solve

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import io

import networkx as nx
import numpy as np
import pygraphviz
from networkx.drawing.nx_agraph import graphviz_layout
import errno, os, signal, functools

def _timeout(timeout_name, error_message=os.strerror(errno.ETIME)):
    '''A timeout decorator that will evaluate class variable for timeout duration'''
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise MOG.TimeoutError(error_message)

        def wrapper(self, *args, **kwargs):
            seconds =  getattr(self, timeout_name)
            if seconds is not None:
                signal.signal(signal.SIGALRM, _handle_timeout)
                signal.alarm(seconds)
            try:
                result = func(self, *args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return functools.wraps(func)(wrapper)

    return decorator
class MOG:

    @property
    def graph(self):
        return self._graph

    @property
    def path_matrix(self):
        return self._path_matrix

    def __init__(self, smiles, symmetry=True, build_timeout=None):
        self._G = smiles2graph(smiles)
        self._graph = nx.DiGraph()
        self._symmetry = symmetry
        self._path_matrix = None
        self._build_timeout = build_timeout
        self._build()

    class TimeoutError(Exception):
        pass

    @_timeout('_build_timeout')
    def _build(self):
        #build the MOG
        #create starting graph, which is quotient graph
        equiv_fxn = lambda a, b: False
        if self._symmetry:
            equiv_fxn = equiv_function(equiv_classes(self._G, 'atom_type', 'bond'))

        qgraph = nx.algorithms.minors.quotient_graph(self._G, equiv_fxn)

        self._label_map = {}
        mapping = {}
        for n in qgraph.nodes():
            if len(n) > 1:
                new_n = set()
                new_n.add(n)
                new_n = frozenset(new_n)
                label = []
                label.append('{{{}}}'.format(','.join(self._G.node[ni]['atom_type'] for ni in n)))
                label.append('{{{}}}'.format(','.join([str(ni) for ni in n])))
                self._label_map[n] = label
                mapping[n] = new_n
            else:
                node = list(n)[0]
                self._label_map[node] = [self._G.node[node]['atom_type'], str(node)]


        nx.relabel_nodes(qgraph, mapping, copy=False)

        #add atoms to MOG
        atom_beads = []
        for n in qgraph:
            self._add_bead(n, self._graph, self._G)
            atom_beads.append(n)
        self._node_layers = [atom_beads]

        #contract edges on quotient graph
        self._node_layers.append(self._remove_bond(qgraph, self._G))

        #agglomerate beads to finish out MOG
        self._agglomerate_layer(self._graph, self._G)
        self._root = self._node_layers[-1][0]
        self._graph.node[self._root]['label'] = 'root'
        self._build_path_matrix()


    def _add_bead(self, nbunch, graph, mol):
        if nbunch is int:
            nbunch = frozenset([nbunch])
        atom_string = ','.join([self._label_map[n][0] for n in nbunch])
        id_string = ','.join([self._label_map[n][1] for n in nbunch])
        if nbunch not in self._graph:
            self._graph.add_node(nbunch, **{'label': atom_string + '\n' + id_string})
            return True
        return False

    def _remove_bond(self, graph, mol):
        #now remove bonds
        new_beads = []
        for bond in graph.edges():
            nbunch = set(bond[0])
            nbunch |= set(bond[1])
            nbunch = frozenset(nbunch)
            if self._add_bead(nbunch, graph, mol):
                self._graph.add_edge(nbunch, bond[0])
                self._graph.add_edge(nbunch, bond[1])
                new_beads.append(nbunch)
        return new_beads

    def __getitem__(self, index):
        return self._node_layers[index]

    def __len__(self):
        return len(self._node_layers)

    def _agglomerate_layer(self, graph, mol):
        nodes = self._node_layers[-1]
        new_nodes = set()
        for i,ni in enumerate(nodes):
            for nj in nodes[i + 1:]:
                if len(ni - nj) == 1 and len(nj - ni) == 1:
                    nbunch = frozenset(ni | nj)
                    new_nodes.add(nbunch)
                    self._add_bead(nbunch, graph, mol)
                    graph.add_edge(nbunch, ni)
                    graph.add_edge(nbunch, nj)
        if len(new_nodes) > 0:
            self._node_layers.append(list(new_nodes))
            self._agglomerate_layer(graph, mol)

    def _build_path_matrix(self):
        '''Build path matrix by exploring all root to child paths'''
        if self._path_matrix is None:
            #reverse index so that root is start of paths
            self._path_map = {k: i for i,k in enumerate(self[0])}
            index = len(self._path_map)
            for n in self._graph.nodes():
                if n not in self._path_map:
                    self._path_map[n] = index
                    index += 1
            #self._path_map = {k: len(self._graph) - i - 1 for i,k in enumerate(self._graph.nodes())}
            self._path_matrix =[[0 for _ in range(len(self._graph))]]
            self._build_path_layer(self._root)
            #make OR operated path_layer
            self._path_valid_matrix = np.zeros( (len(self[0]), len(self._graph)), dtype=np.int32)
            for child in self[0]:
                for path in self._path_matrix:
                    if path[self._path_map[child]] == 1:
                        #in theory, the keys for the atom
                        #nodes should start at 0 so there is no index error
                        self._path_valid_matrix[self._path_map[child]] += path
            #convert from integers to 0/1 indicators
            self._path_valid_matrix = self._path_valid_matrix == 1

    def _build_path_layer(self, node, index = 0):

        self._path_matrix[index][self._path_map[node]] = 1

        if len(self._graph[node]) == 0:
            return
        nlist = list(self._graph.neighbors(node))
        for i,n in enumerate(nlist):
            # if we have more than 1 child, we're creating additional paths
            new_index = index
            #save modifying the original path for last
            if i != len(nlist) - 1:
                #duplicate current path
                self._path_matrix.append(self._path_matrix[index][:])
                new_index = len(self._path_matrix) - 1

            # descend
            self._build_path_layer(n, new_index)

    def mappings(self):
        labels = ['({})'.format(d['label'].split('\n')[0]) for _,d in self._graph.nodes(data=True)]
        for s in integer_solve(self.path_matrix):
            yield [labels[i] for i,si in enumerate(s) if si == 1]

    def print_path(self, p):
        inv_map = {v: k for k, v in self._path_map.items()}
        path_str = ''
        for i,pi in enumerate(p):
            if pi == 1:
                path_str += '(' + self._graph.node[inv_map[i]]['label'].split('\n')[0] + ')-'
        print(path_str[:-1])#omit final dash

    def draw(self, format='svg'):
        pos = graphviz_layout(self._graph, prog='dot', args='')
        fig = plt.figure(figsize=(12, 8))
        if len(self._graph) <= 25:
            nx.draw(self._graph, pos, labels={n: d['label'] for n,d in self._graph.nodes(data=True)},
                node_size=5000)
        else:
            nx.draw(self._graph, pos, index + 1)
        with io.BytesIO() as output:
            fig.savefig(output, format=format)
            fig.clf()
            plt.close('all')
            return output.getvalue()