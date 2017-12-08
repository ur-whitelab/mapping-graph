import networkx as nx
import operator

def equiv_classes(graph, node_key='bond', edge_key='atom_type'):
    '''Function to identify equivalent bonds'''
    #Sub trees are constructed setting each node of edge graph LG to be root
    sub_trees = dict()

    def node_equal(n1, n2):
        '''Determine if two nodes are equal'''
        return n1[node_key] == n2[node_key]
    def edge_equal(e1, e2):
        return e1[edge_key] == e2[edge_key]

    #equivalence classes
    equiv = [set([x]) for x in graph]
    #iterate over the autoisomorphisms
    gm = nx.isomorphism.GraphMatcher(graph, graph,
                                       node_match=node_equal,
                                       edge_match=edge_equal)
    for map in gm.isomorphisms_iter():
        #iterate over the bijective mapping
        for k1,k2 in map.items():
            if k1 == k2:
                continue
            #found a swap of two nodes
            #they must be in the same equivalence class
            for i in range(len(equiv)):
                if k1 in equiv[i]:
                    break
            for j in range(len(equiv)):
                if k2 in equiv[j]:
                    break
            if(i != j):
                equiv[i] |= equiv[j]
                del equiv[j]
    return equiv
