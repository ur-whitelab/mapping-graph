from .mog import MOG
from .reading import smiles2graph, draw, chem_line_graph
from .equiv import equiv_classes
import fire
import networkx as nx
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pygraphviz
from networkx.drawing.nx_agraph import graphviz_layout
import random

import warnings, matplotlib as mpl
warnings.filterwarnings("ignore", category=mpl.mplDeprecation)
warnings.filterwarnings("ignore", category=UserWarning)

def start():
    fire.Fire({
        
        'MOG': mog,
        'draw': draw_mol,
        'count': count
    })

def mog(smiles, output='mog.png', symmetry=True, paths=False):
    mog = MOG(smiles, symmetry)
    #print(len(mog),len(smiles2graph(smiles)))
    if paths:
        for p in mog.path_matrix:
            mog.print_path(p)
    if output is not None:
        plot = mog.draw(format=output.split('.')[1])
        with open(output, 'wb') as f:
            f.write(plot)
def count(smilesdata, output='compression.png',timeout=5,symmetry=True):
    x=[]
    compression=[]
    counter=0
    with open(smilesdata) as infile:
        for line in infile:
            smiles=line.split()[1]
            G= smiles2graph(smiles)
            atom_classes=equiv_classes(G,node_key='atom_type',edge_key='bond')
            #print(G.number_of_edges(),len(atom_classes))
            if (G.number_of_edges()==0) or (len(atom_classes)==1):
                continue
            try:
                mog = MOG(smiles,symmetry,timeout)
            except MOG.TimeoutError:                
                print(smiles)
                continue
            LG=chem_line_graph(G)
            bond_classes = equiv_classes(LG)
            comp_ratio=len(mog.graph)/float(2**len(bond_classes)-1)
            if comp_ratio >=1:
                #print('compression ratio',comp_ratio)
                counter+=1
            x.append(G.number_of_nodes())
            compression.append(comp_ratio)
    x_bins=max(x)-min(x)
    y_bins=(max(compression)-min(compression))*10
    plt.hist2d(x, compression, bins=(x_bins,y_bins),cmap='Blues')
    plt.colorbar()
    #print('counter', counter)
    #plt.scatter(x,compression,c='b',alpha=0.5)        
    plt.savefig(output,format='png')


def draw_mol(smiles, output='molecule.svg', graph=None, line_graph=None):
    G=smiles2graph(smiles)
    LG= chem_line_graph(G)
    bond_classes = equiv_classes(LG)
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
    if graph is not None:
        pos = graphviz_layout(G, prog='neato')
        plt.figure(figsize=(4,4))
        nx.draw(G, pos, node_color=colors, labels={n: d['atom_type'] for n,d in G.nodes(data=True)})
        plt.savefig(graph)

    if line_graph is not None:
        colors = [None for n in LG]
        for i,n in enumerate(LG):
            for j,bc in enumerate(bond_classes):
                if n in bc:
                    colors[i] = cmap(j)
        pos = graphviz_layout(LG, prog='neato')
        plt.figure(figsize=(4,4))
        nx.draw(LG, pos, node_color=colors, labels={n: d['bond'] for n,d in LG.nodes(data=True)})
        plt.savefig(line_graph)


