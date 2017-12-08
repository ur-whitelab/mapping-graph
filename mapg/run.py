from .mot import MOT
from .reading import smiles2graph, draw
from .equiv import equiv_classes
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
        'draw': draw_mol
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
    atom_classes = equiv_classes(G, 'atom_type', 'bond')
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


