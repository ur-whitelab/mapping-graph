def hash_neighs(n, graph, depth=0, tree=None, max_depth=1):
    '''Builds a tree to fixed depth of all neighbors of n'''

    
    if(tree is None):
        tree = nx.Graph(root=graph.node[n]['bond'])
        tree.add_node(n, bond=graph.node[n]['bond'])

    if(depth == max_depth):
        return tree
    for neigh in nx.all_neighbors(graph, n):
        if(not neigh in tree.node):
            tree.add_node(neigh, bond=graph.node[neigh]['bond'])
            tree.add_edge(neigh, n)
    for neigh in nx.all_neighbors(graph, n):
        if(neigh in tree.node):    
        
            tree = hash_neighs(neigh, graph, depth + 1, tree, max_depth)
    return tree

def bond_equiv_classes(G, LG, depth=None):
    '''Function to identify equivalent bonds'''
    if depth==None:
        ''' If depth is not specified, it is set to the 
        height of maximum spanning tree for LG'''
        for g_node in LG.nodes_iter():
            spath=nx.shortest_path_length(LG,source=g_node)
        depth=max(spath.items(), key=operator.itemgetter(1))[1]
        
    #Sub trees are constructed setting each node of edge graph LG to be root
    sub_trees = dict()

    def node_equal(n1, n2):
        '''Determine if two nodes are isomorphic'''
        return (n1['bond'] == n2['bond'])
    
    #build all the trees
    for i,n in enumerate(LG.nodes_iter()): 
        p = hash_neighs(n, LG, max_depth=depth)
        sub_trees[n] = p

    #equivalence classes
    equiv = [set() for x in LG.nodes_iter()]
    for e,n in zip(equiv, LG.nodes_iter()):
        e.add(n)
    
    #find the isomorphic trees and equivalence classes
    
    for k1,g1 in sub_trees.items():
       
        
        for k2, g2 in sub_trees.items():

            if(k1 == k2):
                continue
            if (LG.node[k1]['bond']==LG.node[k2]['bond']):
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
