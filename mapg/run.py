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

def draw_mol(smiles, output='molecule.svg', graph='graph.svg', line_graph='line_graph.svg'):
    G, LG = smiles2graph(smiles)
    bond_classes = equiv_classes(LG)
    print(len(bond_classes))
    svg = draw(smiles, bond_classes, True)
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

    colors = [None for n in LG]
    for i,n in enumerate(LG):
        for j,bc in enumerate(bond_classes):
            if n in bc:
                colors[i] = cmap(j)
    pos = graphviz_layout(LG, prog='neato')
    plt.figure(figsize=(4,4))
    nx.draw(LG, pos, node_color=colors, labels={n: d['bond'] for n,d in LG.nodes(data=True)})
    plt.savefig(line_graph)


def subtrees(smiles, output='subtrees.svg'):
    G, LG = smiles2graph(smiles)
    #print out equiv classes first

    atom_classes = equiv_classes(G, key='atom_type')
    cmap = plt.cm.get_cmap('Accent')
    colors = [None for n in G]
    for i,n in enumerate(G):
        for j,ac in enumerate(atom_classes):
            if n in ac:
                colors[i] = cmap(j)

    #we build a composed tree for layout
    #then plot the individual subtrees
    composed_trees = nx.DiGraph()
    offset_node = 0
    subtrees = []
    for n in G:
        s = hash_neighs(n, G, 'atom_type')
        subtrees.append(s)
        #need to join them for plotting purposes
        nx.relabel_nodes(s, lambda n: n + offset_node, copy=False)
        offset_node += len(s)
        composed_trees.add_nodes_from(composed_trees.nodes(data=True) + s.nodes(data=True))
        composed_trees.add_edges_from(composed_trees.edges() + s.edges())
    fig = plt.figure(1, figsize=(14,8))
    pos = graphviz_layout(composed_trees, prog='dot')
    labels = {n : composed_trees.node[n]['atom_type'] for n in composed_trees}
    # color nodes the same in each connected subgraph
    for g, c, title in zip(subtrees, colors,
                        ['{}-{}'.format(n, d['atom_type']) for n,d in G.nodes(data=True)]):
        nx.draw(g,
             pos,
             node_color=c,
             vmin=0.0,
             vmax=1.0,
             labels=labels)

    plt.savefig(output)
    print(labels)
    for i, e in enumerate(atom_classes):
        print('{}.'.format(i), ','.join([G.node[a]['atom_type'] for a in e]))

