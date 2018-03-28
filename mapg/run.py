from .mog import MOG
from .reading import smiles2graph, draw, chem_line_graph
from .equiv import equiv_classes
import fire
import rdkit
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
plt.switch_backend('agg')
import pygraphviz
from networkx.drawing.nx_agraph import graphviz_layout
import random
from sympy import bell
import warnings, matplotlib as mpl
warnings.filterwarnings("ignore", category=mpl.mplDeprecation)
warnings.filterwarnings("ignore", category=UserWarning)

def start():
    fire.Fire({

        'MOG': mog,
        'draw': draw_mol,
        'count': count
    })

def mog(smiles, output='mog.png', symmetry=True, paths=False, mappings=True):
    '''
    Compute exhaustive mapping graph operator
    '''
    mog = MOG(smiles, symmetry)
    #print(len(mog),len(smiles2graph(smiles)))
    if paths:
        for p in mog.path_matrix:
            mog.print_path(p)
    if output is not None:
        plot = mog.draw(format=output.split('.')[1])
        with open(output, 'wb') as f:
            f.write(plot)
    for mapping in mog.mappings():
        print('-'.join(mapping))

def count(smilesdata, output='lib_data.txt',timeout=5,symmetry=True, sample_size=8000):
    '''
    Count...?
    '''
    x = []
    compression = []
    #counter = 0
    line_count = 0
    #sample_size=int(sample_size)
    for file_lines in open(smilesdata).readlines(): line_count += 1

    if line_count <= sample_size:
        sample_idx = np.arange(line_count)
    else:
        sample_idx=np.random.choice(line_count,sample_size,replace=False)
    f=open(output,"w")
    f.write('#SMILES heavy_atoms bond_number bell_number naive_count starsbars symmetry_count MOG_nodes \n')
    with open(smilesdata) as infile:
        for counter,line in enumerate(infile):
            if counter not in sample_idx:
                continue
            smiles = line.split()[1]
            #print(smiles)
            try:
                G = smiles2graph(smiles)

            except:
                continue
            mol = rdkit.Chem.MolFromSmiles(smiles)
            heavy_atoms = mol.GetNumHeavyAtoms()
            atoms = len(G)
            edges = G.number_of_edges()
            LG = chem_line_graph(G)
            bond_classes = equiv_classes(LG)
            atom_classes = equiv_classes(G,node_key='atom_type',edge_key='bond')
            bellnum = bell(atoms)-1
            naive_count = (2**edges)-1
            product = 1
            for i in range(len(bond_classes)):
                product *= (len(bond_classes[i])+1)
            stars_bars = product-1
            symmetric_counts = (2**(len(bond_classes)))-1
            #print(G.number_of_edges(),len(atom_classes))
            if (edges == 0) or (len(atom_classes) == 1):
                continue
            try:
                mog = MOG(smiles,symmetry,timeout)
                f.write('{smiles} {h_atoms} {bond_number} {bell} {naive} {starsbars} {symmetric} {mog_nodes} \n'.
                        format(smiles = smiles,
                               h_atoms = heavy_atoms,
                               bond_number = edges,
                               bell = bellnum,
                               naive = naive_count,
                               starsbars = stars_bars,
                               symmetric = symmetric_counts,
                               mog_nodes = len(mog.graph)))
            except MOG.TimeoutError:
                f.write('{smiles} {h_atoms} {bond_number} {bell} {naive} {starsbars} {symmetric} {mog_nodes} \n'.
                        format(smiles = smiles,
                               h_atoms = heavy_atoms,
                               bond_number = edges,
                               bell = bellnum,
                               naive = naive_count,
                               starsbars = stars_bars,
                               symmetric = symmetric_counts,
                               mog_nodes = 'NA'))
                continue

    f.close()

def draw_mol(smiles, output='molecule.svg', graph=None, line_graph=None):
    '''
    Draw a molecule, give information about its equivalent atoms and optionally convert to graph/line graph
    '''
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


