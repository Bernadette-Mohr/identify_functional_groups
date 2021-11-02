#
#  Original authors: Richard Hall and Guillaume Godin
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

#
#
# Richard hall 2017
# IFG main code
# Guillaume Godin 2017
# refine output function
# astex_ifg: identify functional groups a la Ertl, J. Cheminform (2017) 9:36
from rdkit import Chem
from collections import namedtuple


def merge(mol, marked, aset):
    bset = set()
    for idx in aset:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            jdx = nbr.GetIdx()
            if jdx in marked:
                marked.remove(jdx)
                bset.add(jdx)
    if not bset:
        return
    merge(mol, marked, bset)
    aset.update(bset)


carbons_dict = dict()
# atoms connected by non-aromatic double or triple bond to any heteroatom
# c=O should not match (see fig1, box 15).  I think using A instead of * should sort that out?
carbons_dict['DOUBLE_TRIPLE'] = Chem.MolFromSmarts('A=,#[!#6]')
# atoms in non aromatic carbon-carbon double or triple bonds
carbons_dict['CC_DOUBLE_TRIPLE'] = Chem.MolFromSmarts('C=,#C')
# acetal carbons, i.e. sp3 carbons connected to tow or more oxygens, nitrogens or sulfurs; these O, N or S atoms must
# have only single bonds
carbons_dict['ACETAL'] = Chem.MolFromSmarts('[CX4](-[O,N,S])-[O,N,S]')
# all atoms in oxirane, aziridine and thiirane rings
carbons_dict['OXIRANE_ETC'] = Chem.MolFromSmarts('[O,N,S]1CC1')


def get_pattern_name(patterns, group):
    for key, value in patterns.items():
        if all([v in group for v in value]):
            return key


def identify_functional_groups(mol):

    patterns = dict()
    marked = set()
    # mark all heteroatoms in a molecule, including halogens
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 1):  # would we ever have hydrogen?
            marked.add(atom.GetIdx())

    # mark the four specific types of carbon atom
    for name, patt in carbons_dict.items():
        if mol.GetSubstructMatch(patt):
            for path in mol.GetSubstructMatches(patt):
                atomindices = list()
                for atomindex in path:
                    atomindices.append(atomindex)
            # not ideal, multiple occurrences of the same pattern aren't captured with the dictionary!
            patterns[name] = atomindices
            marked.update(atomindices)

    # merge all connected marked atoms to a single FG
    groups = list()
    while marked:
        grp = {marked.pop()}
        merge(mol, marked, grp)
        groups.append(grp)

    # extract also connected unmarked carbon atoms
    ifg = namedtuple('IFG', ['atomIds', 'atoms', 'type', 'name'])
    ifgs = []
    for g in groups:
        # print(g)
        uca = set()
        for atomidx in g:
            for n in mol.GetAtomWithIdx(atomidx).GetNeighbors():
                if n.GetAtomicNum() == 6:
                    uca.add(n.GetIdx())
        ifgs.append(ifg(atomIds=tuple(list(g)), atoms=Chem.MolFragmentToSmiles(mol, g, canonical=True),
                        type=Chem.MolFragmentToSmiles(mol, g.union(uca), canonical=True),
                        name=get_pattern_name(patterns, g)))
    return ifgs


if __name__ == "__main__":
    import sys

    identify_functional_groups(Chem.MolFromSmiles(sys.argv[1]))
