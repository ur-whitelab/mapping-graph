import networkx as nx

def merge_nodes(G, nodes, new_node, fg2cg=dict(), attr_dict=None, **attr):
    """
    Merges the selected `nodes` of the graph G into one `new_node`,
    meaning that all the edges that pointed to or from one of these
    `nodes` will point to or from the `new_node`.
    attr_dict and **attr are defined as in `G.add_node`.
    """
    #print('merging {} into {} on graph nodes {} \n with map {}'.format(nodes, new_node, G.nodes(), fg2cg))
    G.add_node(new_node, attr_dict, **attr) # Add the 'merged' node

    for n1,n2,data in G.edges(data=True):
        # For all edges related to one of the nodes to merge,
        # make an edge going to or coming from the `new gene`.
        #convert them
        if n1 in nodes:
            G.add_edge(new_node,n2,data)
        elif n2 in nodes:
            G.add_edge(n1,new_node,data)

    for n in nodes: # remove the merged nodes
        try:
            G.remove_node(n)
        except nx.NetworkXError:
            print('error while deleting {}, here are nodes {}'.format(n, G.nodes()))
            print('FG2CG', fg2cg)


    for n in nodes:
        fg2cg[n] = new_node


def expand_from_edgeset(fine_graph, deleted_edges):
    #print('-----Building new graph------')
    G = fine_graph.copy()
    #add size attribute
    for n,d in G.nodes(data=True):
        d['size'] = 1
    fg2cg = dict()
    for de in deleted_edges:

        #need to convert delete_edge node references into current graph node-set
        #while loops are because we can have 3 merging with 4 into CG-0 which
        #merges with CG-1 which merges with CG-4. So if we want find 3
        #in current graph, we have to lookup 3, then lookup CG-0 then look up CG-1 etc
        nodes = []
        for e in de:
            while e in fg2cg:
                e = fg2cg[e]
            nodes.append(e)
        #make new node's dictionary
        attr = dict()
        attr['size'] = sum([G.node[n]['size'] for n in nodes])
        attr['atmType'] = ','.join([G.node[n]['atmType'] for n in nodes])
        merge_nodes(G, nodes, 'CG-{}'.format(fine_graph.number_of_nodes() - G.number_of_nodes()), fg2cg=fg2cg, attr_dict=attr)
    return G, fg2cg