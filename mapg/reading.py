'''Functions to help read in molecular structures'''
from rdkit import Chem
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pygraphviz
import matplotlib as mpl
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
import svgutils.transform as sg


def smiles2graph(sml):
    '''Argument for the RD2NX function should be a valid SMILES sequence'''
    m = Chem.MolFromSmiles(sml)
    m = Chem.AddHs(m)
    G = nx.Graph()

    for i in m.GetAtoms():
        G.add_node(i.GetIdx(),atmType=i.GetSymbol())

    for j in m.GetBonds():
        G.add_edge(j.GetBeginAtomIdx(),j.GetEndAtomIdx())

    #get line graph (vertx -> edge, edge -> vertex)
    LG = nx.line_graph(G)
    #add the data to edges for atom types. Note, doesn't include bond order
    for n in LG.nodes():
        LG.node[n]['bond'] = [G.node[n[0]]['atmType'], G.node[n[1]]['atmType']]
        LG.node[n]['bond'].sort()
        LG.node[n]['bond'] = ''.join(LG.node[n]['bond'])
    return G, LG

def draw(sml, equiv_bonds=None, color_by_equiv=False):
    '''Draw a structure with equivalent bonds optionally highlighted from a SMILES string'''
    m = Chem.MolFromSmiles(sml)
    m = Chem.AddHs(m)
    rdDepictor.Compute2DCoords(m)


    drawer = rdMolDraw2D.MolDraw2DSVG(400,200)
    drawer.drawOptions().bgColor = None
    if equiv_bonds is not None:
        #convert to list
        def lookup_bond(b):
            for i,e in enumerate(equiv_bonds):
                #print(b.GetBeginAtomIdx(), b.GetEndAtomIdx(),e)
                if( (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) in e):

                    #only include bond classes that have more than one member
                    if(len(e) > 1):
                        return i

                    else:
                        break
            return -1
        bond_classes = [lookup_bond(b) for b in m.GetBonds()]
        #print(bond_classes)
        #remove the -1s
        highlight = []
        highlight_atoms = set()
        for i, b in enumerate(bond_classes):
            if b >= 0:
                highlight.append(i)
                highlight_atoms.add(m.GetBonds()[i].GetBeginAtomIdx())
                highlight_atoms.add(m.GetBonds()[i].GetEndAtomIdx())
        bn = len(bond_classes)
        ##cmap=plt.cm.get_cmap('Accent', bn)##
        cmap=plt.cm.get_cmap('Accent')
        #print(bn)
        #use none or cmap for the classes
        colors={ i : cmap(bond_classes[i]) for i in highlight}
        ##colors={ i : cmap(bond_classes[i] / bn) for i in highlight}##
        #print(colors)
        if(color_by_equiv):
            drawer.DrawMolecule(m,highlightAtoms=list(highlight_atoms),
                            highlightBonds=highlight,
                            highlightBondColors=colors,
                            highlightAtomColors={k: (0.8,0.8,0.8) for k in list(highlight_atoms)})

        else:
             drawer.DrawMolecule(m,highlightAtoms=list(highlight_atoms),
                            highlightBonds=highlight)
    else:
        drawer.DrawMolecule(m)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    if os.path.isdir("./img") == False:
        os.mkdir('./img')
    target = open('./img/SVG'+sml+'.svg', 'w')
    target.write(svg)
    target.close()
    return display.SVG(svg)

