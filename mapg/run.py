from .mot import MOT
from .reading import smiles2graph, draw
from .equiv import equiv_classes, hash_neighs
import fire
import networkx as nx
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pygraphviz
from networkx.drawing.nx_agraph import graphviz_layout
import random

def start():
    fire.Fire({
        'mot': mot,
        'draw': draw_mol,
        'subtrees': subtrees
    })

def mot(smiles, symmetry=True,mot_output='mot.svg'):
    mot = MOT(smiles, symmetry)
    mot.build()
    mot.prune_parents()
    mot.prune_nodes()
    if mot_output is not None:
        svg = mot.draw()
        with open(mot_output, 'wb') as f:
            f.write(svg)

def draw_mol(smiles, output='molecule.svg', graph='graph.svg'):
    G, LG = smiles2graph(smiles)
    bec = equiv_classes(LG)
    svg = draw(smiles, bec, True)
    with open(output, 'w') as f:
        f.write(svg)
    atom_classes = equiv_classes(G, key='atom_type')
    cmap = plt.cm.get_cmap('Accent')
    colors = [None for n in G]
    for i,n in enumerate(G):
        for j,ac in enumerate(atom_classes):
            if n in ac:
                colors[i] = cmap(j)
    pos = graphviz_layout(G, prog='neato')
    plt.figure(figsize=(4,4))
    nx.draw(G, pos, node_color=colors, labels={n: d['atom_type'] for n,d in G.nodes(data=True)})
    plt.savefig(graph)


def subtrees(smiles, output='subtrees.svg'):
    from collections import deque
    G, LG = smiles2graph(smiles)
    #we build a composed tree for layout
    #then plot the individual subtrees
    composed_trees = nx.DiGraph()
    offset_node = 0
    subtrees = []
    for n in G:
        s = hash_neighs(deque([n]), G, 'atom_type')
        subtrees.append(s)
        #need to join them for plotting purposes
        nx.relabel_nodes(s, lambda n: n + offset_node, copy=False)
        offset_node += len(s)
        composed_trees.add_nodes_from(composed_trees.nodes(data=True) + s.nodes(data=True))
        composed_trees.add_edges_from(composed_trees.edges() + s.edges())
    fig = plt.figure(1, figsize=(8,8))
    pos = graphviz_layout(composed_trees, prog='dot')
    labels = {n : composed_trees.node[n]['atom_type'] for n in composed_trees}
    # color nodes the same in each connected subgraph
    for g, title in zip(subtrees,
                        ['{}-{}'.format(n, d['atom_type']) for n,d in G.nodes(data=True)]):
        c = [ random.random() ] * len(g) # random color...
        nx.draw(g,
             pos,
             node_color=c,
             vmin=0.0,
             vmax=1.0,
             labels=labels)

    plt.savefig(output)
    for i, e in enumerate(equiv_classes(G, key='atom_type')):
        print('{}.'.format(i), ','.join([labels[a] for a in e]))
