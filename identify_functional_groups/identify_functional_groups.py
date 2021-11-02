import pandas as pd
import argparse
import configparser
from pathlib import Path
from rdkit import Chem
from ifg import identify_functional_groups
from elements import get_element_names

# configparser object made available globally, because passing as argument in pandas.apply() not possible.
fg_patterns = configparser.ConfigParser()
fg_names = configparser.ConfigParser()


def get_fgs(pattern_file, names_file):
    fg_patterns.read(pattern_file)
    fg_names.read(names_file)


def isRingAromatic(mol, ring):
    for idx in ring:
        if not mol.GetBondWithIdx(idx).GetIsAromatic():
            return False
        else:
            return True


def getRingSystems(mol, includeSpiro=False):
    rings = mol.GetRingInfo()
    systems = []
    aromaticity = []
    for a_ring, b_ring in zip(rings.AtomRings(), rings.BondRings()):
        ringAts = set(a_ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and (includeSpiro or nInCommon > 1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        aromatic = isRingAromatic(mol, b_ring)
        # print(aromatic)
        aromaticity.append(aromatic)
        nSystems.append(ringAts)
        systems = nSystems
    # print(aromaticity)
    if len(systems) == len(aromaticity):
        # print(f'IF')
        systems = list(zip(systems, aromaticity))
    else:
        # print(f'ELSE')
        if any(ar is True for ar in aromaticity):
            # print('True')
            systems = list(zip(systems, [True] * len(systems)))
        else:
            # print('False')
            systems = list(zip(systems, [False] * len(systems)))
    # print(systems)
    return systems


def select_most_important(recognized_fgs):
    fallback_groups = [['not identified'], ['halogen deriv.'], ['alkene', 'alkyne'],
                       ['carbonyl group', 'methoxy', 'ethoxy', 'carboxylic acid deriv.'],
                       ['carbonyl with nitrogen', 'amine', 'hydrazine deriv.', 'carbamic acid deriv.'],
                       ['4-atom-heterocycle'], ['5-atom-heterocycle'], ['6-atom-heterocycle']]
    # print('called!')
    # print(recognized_fgs)
    if len(recognized_fgs) > 1:
        for group in fallback_groups:
            while any([fg in recognized_fgs for fg in group]) and len(recognized_fgs) > 1:
                for fg in group:
                    if len(recognized_fgs) > 1:
                        recognized_fgs.discard(fg)

    return recognized_fgs


def select_important_rings(recognized_fgs):
    fallback_groups = [['4-atom-heterocycle'], ['5-atom-heterocycle'], ['6-atom-heterocycle']]
    # print('called!')
    # print(recognized_fgs)
    if len(recognized_fgs) > 1:
        for group in fallback_groups:
            while any([fg in recognized_fgs for fg in group]) and len(recognized_fgs) > 1:
                for fg in group:
                    if len(recognized_fgs) > 1:
                        recognized_fgs.discard(fg)

    return recognized_fgs


def exclude_ring_elements(mol, pattern_name, keys, ring=None):
    if ring is not None:
        recognized_fgs = list()
        for key in keys:
            pattern = fg_patterns.get(pattern_name, key)
            for path in mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)):
                if not all([atomId in ring[0] for atomId in path]):
                    recognized_fgs.append(fg_names.get(pattern_name, key))
        return recognized_fgs

    else:
        return [fg_names.get(pattern_name, k) for k in keys
                if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get(pattern_name, k)))]


def process_multiple_atoms(mol, hetatms, elem_names, ring=None):
    recognized_fgs = set()
    # print('inside process_multiple_atoms:', hetatms)
    section_name = '_AND_'.join([elem_names[hetatm] for hetatm in hetatms])
    # print(section_name)
    keys = fg_patterns.options(section_name)
    # print([fg_names.get(section_name, k) for k in keys
    #        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get(section_name, k)))])
    if ring:
        recognized_fgs.update(exclude_ring_elements(mol, section_name, keys, ring))
    else:
        recognized_fgs.update(exclude_ring_elements(mol, section_name, keys))

    return recognized_fgs


def process_fragments(mol, fg, elem_names, recognized_fgs, ring=None):
    hetatms = set()
    hetatms.update([char.upper() for char in fg.atoms if char.isalpha()
                    and char.casefold() != 'c' and char.casefold() != 'h'])
    hetatms = sorted(hetatms)

    if len(hetatms) == 0:
        if 'CC_DOUBLE_TRIPLE' in fg.name:
            pattern_name = 'CC_DOUBLE_TRIPLE'
            keys = fg_patterns.options(pattern_name)
            # print([fg_names.get('CC_DOUBLE_TRIPLE', k) for k in keys
            #        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get('CC_DOUBLE_TRIPLE', k)))])
            if ring:
                recognized_fgs.update(exclude_ring_elements(mol, pattern_name, keys, ring))
            else:
                recognized_fgs.update(exclude_ring_elements(mol, pattern_name, keys))
        else:
            pattern_name = 'ALKANE'
            keys = fg_patterns.options(pattern_name)
            # print([fg_names.get('ALKANE', k) for k in keys
            #        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get('ALKANE', k)))])
            if ring:
                recognized_fgs.update(exclude_ring_elements(mol, pattern_name, keys, ring))
            else:
                recognized_fgs.update(exclude_ring_elements(mol, pattern_name, keys))
    else:
        if any(elem in hetatms for elem in ['F', 'Cl', 'Br', 'I']):
            pattern_name = 'HALOGEN'
            keys = fg_patterns.options(pattern_name)
            # print([fg_names.get('HALOGEN', k) for k in keys
            #        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get('HALOGEN', k)))])
            if ring:
                recognized_fgs.update(exclude_ring_elements(mol, pattern_name, keys, ring))
            else:
                recognized_fgs.update(exclude_ring_elements(mol, pattern_name, keys))
        elif fg.name and 'ACETAL' in fg.name:
            pattern_name = 'ACETAL'
            keys = fg_patterns.options(pattern_name)
            # print([fg_names.get('ACETAL', k) for k in keys
            #        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get('ACETAL', k)))])
            if ring:
                recognized_fgs.update(exclude_ring_elements(mol, pattern_name, keys, ring))
            else:
                recognized_fgs.update(exclude_ring_elements(mol, pattern_name, keys))
            # all_recognized_fgs.update(recognized_fgs)
        else:
            if len(hetatms) > 1:
                while len(recognized_fgs) == 0:
                    if ring:
                        recognized_fgs.update(process_multiple_atoms(mol, hetatms, elem_names, ring))
                        if len(recognized_fgs) == 0:
                            for atom in hetatms:
                                hetatms.remove(atom)
                                if len(hetatms) != 0:
                                    element = elem_names[atom]
                                    keys = fg_patterns.options(element)
                                    # print([fg_names.get(element, k) for k in keys
                                    #        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get(element, k)))])
                                    recognized_fgs.update(exclude_ring_elements(mol, element, keys, ring))
                                    if len(hetatms) > 1:
                                        recognized_fgs.update(process_multiple_atoms(mol, hetatms, elem_names, ring))
                                else:
                                    recognized_fgs.add('not identified')
                    else:
                        recognized_fgs.update(process_multiple_atoms(mol, hetatms, elem_names))
                        if len(recognized_fgs) == 0:
                            for atom in hetatms:
                                hetatms.remove(atom)
                                if len(hetatms) != 0:
                                    element = elem_names[atom]
                                    keys = fg_patterns.options(element)
                                    # print([fg_names.get(element, k) for k in keys
                                    #        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get(element, k)))])
                                    recognized_fgs.update(exclude_ring_elements(mol, element, keys))
                                    if len(hetatms) > 1:
                                        recognized_fgs.update(process_multiple_atoms(mol, hetatms, elem_names))
                                else:
                                    recognized_fgs.add('not identified')
            else:
                element = elem_names[hetatms.pop()]
                keys = fg_patterns.options(element)
                # print([fg_names.get(element, k) for k in keys
                #        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get(element, k)))])
                # recognized_fgs.update([fg_names.get(element, k) for k in keys
                #                        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get(element, k)))])
                if ring:
                    recognized_fgs.update(exclude_ring_elements(mol, element, keys, ring))
                else:
                    recognized_fgs.update(exclude_ring_elements(mol, element, keys))
    # print(recognized_fgs)
    recognized_fgs = select_most_important(recognized_fgs)
    # print(recognized_fgs)
    return recognized_fgs


def identify_rings(mol, rings, elem_names, fg=None):
    recognized_fgs = set()
    for ring in rings:
        if fg:
            if not all([atomId in _ring[0] for atomId in fg.atomIds for _ring in rings]):
                recognized_fgs = process_fragments(mol, fg, elem_names, recognized_fgs, ring)
                # print(recognized_fgs)
                recognized_fgs = select_most_important(recognized_fgs)
                # print(recognized_fgs)
            if (fg.name and 'OXIRANE_ETC' in fg.name) \
                    or (len(ring[0]) == 3 and any([mol.GetAtomWithIdx(idx).GetSymbol() != 'C' for idx in ring[0]])):
                keys = fg_patterns.options('OXIRANE_ETC')
                # print([fg_names.get('OXIRANE_ETC', k) for k in keys
                #        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get('OXIRANE_ETC', k)))])
                recognized_fgs.update([fg_names.get('OXIRANE_ETC', k) for k in keys
                                       if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get('OXIRANE_ETC', k)))])
            # else:
            size = len(ring[0])
            # print('ring size:', size)
            if size < 7:
                section_name = f'{str(size)}_ATOMS_CYCLE'
                keys = fg_patterns.options(section_name)
                # int([fg_names.get(section_name, k) for k in keys
                #      if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get(section_name, k)))])
                recognized_fgs.update([fg_names.get(section_name, k) for k in keys
                                       if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get(section_name, k)))])
                # recognized_fgs = select_important_rings(recognized_fgs)
            # else:
                # section_name = f'MULTI_RING'
                # keys = fg_patterns.options(section_name)
                # # int([fg_names.get(section_name, k) for k in keys
                # #      if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get(section_name, k)))])
                # recognized_fgs.update([fg_names.get(section_name, k) for k in keys
                #                        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get(section_name, k)))])
                # # recognized_fgs = select_important_rings(recognized_fgs)
        else:
            size = len(ring[0])
            # print('ring size:', size)
            if size < 7:
                section_name = f'{str(size)}_ATOMS_CYCLE'
                # print(section_name)
                keys = fg_patterns.options(section_name)
                # print([fg_names.get(section_name, k) for k in keys
                #        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get(section_name, k)))])
                recognized_fgs.update([fg_names.get(section_name, k) for k in keys
                                       if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get(section_name, k)))])
                # recognized_fgs = select_important_rings(recognized_fgs)

    return recognized_fgs


def store_results(all_recognized_fgs, recognized_fgs):

    if recognized_fgs:
        if all([rfg not in fgs for rfg in recognized_fgs for fgs in all_recognized_fgs]):
            # if all([recognized_fgs.isdisjoint(fgs) for fgs in all_recognized_fgs]):
            for elem in list(recognized_fgs):
                all_recognized_fgs.append(elem)

    return all_recognized_fgs


def identify_patterns(mol, smiles, rings, functional_groups):
    print(f'\n{smiles}')

    elem_names = get_element_names()
    all_recognized_fgs = list()

    if not functional_groups:
        if rings:
            recognized_fgs = identify_rings(mol, rings, elem_names)
        else:
            recognized_fgs = set()
            keys = fg_patterns.options('ALKANE')
            # print([fg_names.get('ALKANE', k) for k in keys
            #        if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get('ALKANE', k)))])
            recognized_fgs.update([fg_names.get('ALKANE', k) for k in keys
                                   if mol.GetSubstructMatch(Chem.MolFromSmarts(fg_patterns.get('ALKANE', k)))])
        all_recognized_fgs = store_results(all_recognized_fgs, recognized_fgs)
    else:
        for fg in functional_groups:
            # print(fg)
            recognized_fgs = set()
            if rings:
                recognized_fgs = identify_rings(mol, rings, elem_names, fg)
                all_recognized_fgs = store_results(all_recognized_fgs, recognized_fgs)
            else:
                recognized_fgs = process_fragments(mol, fg, elem_names, recognized_fgs)
                all_recognized_fgs = store_results(all_recognized_fgs, recognized_fgs)

    if not all_recognized_fgs:
        all_recognized_fgs.append({'not identified'})
    print(f'end result: {all_recognized_fgs}\n')

    return all_recognized_fgs


def analyze_fragments(smiles):
    s = pd.Series(dtype='object')
    mol = Chem.MolFromSmiles(smiles)

    rings = getRingSystems(mol)
    if rings:
        if any(filter(lambda x: x[1] is True, rings)):
            s['aromatic'] = True
        else:
            s['aromatic'] = False

    fg = identify_functional_groups(mol)
    # print(fg)
    s['fg_names'] = identify_patterns(mol, smiles, rings, fg)
    s['IFG'] = [group._asdict() for group in fg]

    return s


def main(path, dataframes, patterns, names):
    get_fgs(patterns, names)

    for frame in dataframes:
        df = pd.read_pickle(frame)
        fgs_df = pd.DataFrame(columns=['SMILES', '5T_CG'])
        if 'CG' in df.columns:
            fgs_df[['SMILES', '5T_CG']] = df[['SMILES', 'CG']]
        else:
            fgs_df = df
        # fgs_df[['n_rings', 'aromatic']]
        # fgs_df = pd.concat([fgs_df, df['SMILES'].swifter.allow_dask_on_strings(enable=True).apply(analyze_fragments)],
        #                    axis=1)
        fgs_df = pd.concat([fgs_df, df['SMILES'].apply(analyze_fragments)], axis=1)
        save_path = path / frame.name
        # print(fgs_df)
        print(save_path)
        fgs_df.to_pickle(save_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Automated recognition of functional groups from SMILES strings.')
    parser.add_argument('-dir', '--directory', type=Path, required=True, help='Directory the results will be saved in.')
    parser.add_argument('-df', '--dataframes', type=Path, nargs='+', required=True,
                        help='Dataframes with molecule info')
    parser.add_argument('-pts', '--patterns', type=Path, help='File with smarts patterns to identify functional groups',
                        default=Path('/media/bernadette/ElementsSE/Homeoffice/PycharmProjects/ChemSpaceMapping/FunctionalGroups.cnf'))
    parser.add_argument('-n', '--names', type=Path, help='File with strings to name functional groups',
                        default=Path('/media/bernadette/ElementsSE/Homeoffice/PycharmProjects/ChemSpaceMapping/FG_names.cnf'))

    args = parser.parse_args()
    main(args.directory, args.dataframes, args.patterns, args.names)
