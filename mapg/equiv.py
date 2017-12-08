import networkx as nx
import operator

def hash_neighs(node, graph, trait):
    '''This wrapper choses two possible search orders, which
    is necessary for treating cycles'''
    g1 = _hash_neighs([node], graph, trait, True)
    g2 = _hash_neighs([node], graph, trait, False)
    h1 = ''.join([graph.node[x][trait] for x in nx.dag_longest_path(g1)])
    h2 = ''.join([graph.node[x][trait] for x in nx.dag_longest_path(g2)])
    if h1 < h2:
        return g1
    return g2

def _hash_neighs(queue, graph, trait, reverse_order=False, tree=None):
    '''Builds a tree to fixed depth of all neighbors of n'''
    n = queue.pop()

    if(tree is None):
        tree = nx.DiGraph(root=graph.node[n][trait])
        tree.add_node(n, attr_dict={trait: graph.node[n][trait]})

    for neigh in sorted(nx.all_neighbors(graph, n), reverse=reverse_order):
        if(not neigh in tree.node):
            tree.add_node(neigh, attr_dict={trait: graph.node[neigh][trait]})
            tree.add_edge(n, neigh)
            queue.append(neigh)

    if(len(queue) == 0):
        return tree
    else:
        return(_hash_neighs(queue, graph, trait, reverse_order, tree))

def equiv_classes(graph, key='bond'):
    '''Function to identify equivalent bonds'''
    #Sub trees are constructed setting each node of edge graph LG to be root
    sub_trees = dict()

    def node_equal(n1, n2):
        '''Determine if two nodes are equal'''
        return (n1[key] == n2[key])

    #build all the trees
    for n in graph:
        p = hash_neighs(n, graph, key)
        sub_trees[n] = p

    #equivalence classes
    equiv = [set([x]) for x in graph]
    #find the isomorphic trees and equivalence classes
    for k1, g1 in sub_trees.items():
        for k2, g2 in sub_trees.items():
            if(k1 == k2):
                break
            if (graph.node[k1][key] == graph.node[k2][key]):
                '''If the root labels are different
                the sub-trees are not isomorphic'''
                gm = nx.isomorphism.DiGraphMatcher(g1, g2,
                                                   node_match=node_equal)
                if(gm.is_isomorphic()):
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
