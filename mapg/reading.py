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
    order_string = {rdkit.Chem.rdchem.BondType.SINGLE: '-',
                    rdkit.Chem.rdchem.BondType.DOUBLE: '=',
                    rdkit.Chem.rdchem.BondType.TRIPLE: '#',
                    rdkit.Chem.rdchem.BondType.AROMATIC: ':'}

    for i in m.GetAtoms():
        G.add_node(i.GetIdx(),atom_type=i.GetSymbol())

    for j in m.GetBonds():
        u = min(j.GetBeginAtomIdx(),j.GetEndAtomIdx())
        v = max(j.GetBeginAtomIdx(),j.GetEndAtomIdx())
        order = j.GetBondType()
        if order in order_string:
            order = order_string[order]
        G.add_edge(u, v, bond='{}{}{}'.format(G.node[u]['atom_type'],order,G.node[v]['atom_type']))
    return G

def chem_line_graph(graph):
    #get line graph (vertx -> edge, edge -> vertex)
    LG = nx.line_graph(graph)
    #add the data to edges for atom types
    for n in LG.nodes():
        LG.node[n]['bond'] = graph[n[0]][n[1]]['bond']
    for e in LG.edges():
        n = e[0][0]
        LG[e[0]][e[1]]['atom_type'] = graph.node[n]['atom_type']
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
                bond_tup=tuple(sorted([b.GetBeginAtomIdx(), b.GetEndAtomIdx()]))
                if( bond_tup in e):
                    #only include bond classes that have more than one member
                    if(len(e) > 1):
                        return i
                    else:
                        break
            return -1
        bond_classes = [lookup_bond(b) for b in m.GetBonds()]
        #remove the -1s
        highlight = []
        highlight_atoms = set()
        for i, b in enumerate(bond_classes):
            if b >= 0:
                highlight.append(i)
                highlight_atoms.add(m.GetBonds()[i].GetBeginAtomIdx())
                highlight_atoms.add(m.GetBonds()[i].GetEndAtomIdx())
        cmap=plt.cm.get_cmap('Accent')
        #use none or cmap for the classes
        colors={ i : cmap(bond_classes[i]) for i in highlight}
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

