'''Functions to help read in molecular structures'''
import rdkit, rdkit.Chem, rdkit.Chem.rdDepictor, rdkit.Chem.Draw
import networkx as nx
import matplotlib.pyplot as plt
import svgutils.transform as sg


def smiles2graph(sml):
    '''Argument for the RD2NX function should be a valid SMILES sequence

    returns: the graph and it's linegraph
    '''
    m = rdkit.Chem.MolFromSmiles(sml)
    m = rdkit.Chem.AddHs(m)
    G = nx.Graph()

    for i in m.GetAtoms():
        G.add_node(i.GetIdx(),atom_type=i.GetSymbol())

    for j in m.GetBonds():
        G.add_edge(j.GetBeginAtomIdx(),j.GetEndAtomIdx())
    return G, chem_line_graph(G)

def chem_line_graph(graph):
    #get line graph (vertx -> edge, edge -> vertex)
    LG = nx.line_graph(graph)
    #add the data to edges for atom types. Note, doesn't include bond order
    for n in LG.nodes():
        LG.node[n]['bond'] = [graph.node[n[0]]['atom_type'], graph.node[n[1]]['atom_type']]
        LG.node[n]['bond'].sort()
        LG.node[n]['bond'] = ''.join(LG.node[n]['bond'])
    return LG

def draw(sml, equiv_bonds=None, color_by_equiv=False):
    '''Draw a structure with equivalent bonds optionally highlighted from a SMILES string'''
    m = rdkit.Chem.MolFromSmiles(sml)
    m = rdkit.Chem.AddHs(m)

    rdkit.Chem.rdDepictor.Compute2DCoords(m)


    drawer = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG(400,200)
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
    return svg

