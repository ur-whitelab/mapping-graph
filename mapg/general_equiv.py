import networkx as nx
import operator
from collections import deque
def general_hash_neighs(queue, graph,trait,tree=None):
    '''Builds a tree to fixed depth of all neighbors of n'''
    n=queue.popleft()

    if(tree is None):

        tree = nx.Graph(root=graph.node[n][trait])
        tree.add_node(n, trait=graph.node[n][trait])

    for neigh in sorted(nx.all_neighbors(graph, n)):
        if(not neigh in tree.node):
            
            tree.add_node(neigh, trait=graph.node[neigh][trait])
            tree.add_edge(neigh, n)
            queue.append(neigh)
    
    if(len(queue) == 0):
        return tree
    else:
        return(general_hash_neighs(queue,graph,trait,tree))   

def general_equiv_classes(G, LG,key='bond'):
    '''Function to identify equivalent bonds'''
    if key=='bond':
        graph=LG
        trait='bond'
    elif key=='atom':
        graph=G
        print(graph.nodes(data=True))
        trait='atom_type'
    else:
        print('Invalid key-type')
    #Sub trees are constructed setting each node of edge graph LG to be root
    sub_trees = dict()

    def node_equal(n1, n2):
        '''Determine if two nodes are isomorphic'''
        return (n1['trait'] == n2['trait'])

    #build all the trees
    for i,n in enumerate(graph.nodes_iter()):
        p = general_hash_neighs(deque([n]), graph,trait)
        sub_trees[n] = p

    #equivalence classes
    equiv = [set() for x in graph.nodes_iter()]
    for e,n in zip(equiv, graph.nodes_iter()):
        e.add(n)

    #find the isomorphic trees and equivalence classes

    for k1,g1 in sub_trees.items():


        for k2, g2 in sub_trees.items():

            if(k1 == k2):
                continue
            if (graph.node[k1][trait]==graph.node[k2][trait]):
                '''If the root labels are different
                the sub-trees are not isomorphic'''
                gm = nx.isomorphism.GraphMatcher(g1, g2, node_match=node_equal)
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
