from collections import defaultdict
import numpy as np
from math import e, pi, sqrt
from scipy.integrate import quad
import gzip
import io
import logging
from tqdm import tqdm
import pandas as pd

from prody import parsePDB
from prody.atomic.hierview import HierView
from prody.atomic import AAMAP

from chimera.utils import vdw_radius

# Constant factor used for PDF calculation for a normal distribution. Pre-computed here once.
NORM_PDF_FACTOR = 1.0/sqrt(2*pi)

tol = np.finfo(float).eps ** 0.25
logger = logging.getLogger(__name__)

# TODO: Whether to ignore residue insertion code (col 27 in pdb) when identifying unique residues
# Should ideally be False, but set to True for backwards compatibility
IGNORE_INSERTION_CODE = True

# TODO: Whether to 'cache' normal distribution overlap numbers, but in reverse order, emulating
# behavior of old code
EMULATE_CACHING_BUG = True

ANNOTATION_COLUMNS = [
    'pdb_id',
    'pdb_chain',
    'resolution',
    'binding_site_id',
    'ligand_id',
    'ligand_chain',
    'ligand_serial_number',
    'binding_site_residue__pdb_index',
    'binding_site_residue__one_index',
    'catalytic_site_residue__pdb_index',
    'catalytic_site_residue__one_index',
    'enzyme_commission_number',
    'go_terms',
    'binding_affinity__manual',
    'binding_affinity__binding_moad',
    'binding_affinity__pdb_bind_cn',
    'binding_affinity__binding_db',
    'uniprot_id',
    'pubmed_id',
    'full_receptor_sequence'
]


def distance_df_to_csv(df, pdb_id, file_path, compressed=False):

    def _format_integral(x):
        s = f'{x:.7g}'
        if s == '0.0001':
            return '1e-04'
        else:
            return s

    def _format_error(x):
        s = f'{x:.2g}'
        if s == '0.0001':
            return '1e-04'
        else:
            return s

    _open = open
    if compressed:
        _open = gzip.open
        if not file_path.endswith('.gz'):
            file_path += '.gz'

    f = _open(file_path, 'wb')
    f.write(
        ('\n'.join(['# All pairwise distances between receptor protein chain residue atoms ' +
                    'and ligand atoms for ' + pdb_id,
                    '# NOTE: columns 8-11 contain the overlap area between Gaussian ' +
                    'distributions centered at each atom with ',
                    '#   standard deviations set to either the van der Waals radii of ' +
                    'the two atoms or to 1.5']) + '\n').encode('utf8'))
    f.write(('\t'.join(['#pdbID-pdbChain', 'receptor_aa_1-index', 'receptor_aa_value',
                        'receptor_atom_id', 'receptor_atom_value', 'ligand_id',
                        'ligand_atom_value',
                        'euclidean_distance',
                        'overlap_vdw_radii', 'integral_error_vdw_radii',
                        'overlap_1.5', 'integral_error_1.5',
                        'full_receptor_sequence']) + '\n').encode('utf8'))

    df = df.copy()
    df['overlap_vdw_radii'] = df.apply(lambda row: _format_integral(row['overlap_vdw_radii']), axis=1)
    df['integral_error_vdw_radii'] = df.apply(lambda row: _format_error(row['integral_error_vdw_radii']), axis=1)
    df['overlap_1.5'] = df.apply(lambda row: _format_integral(row['overlap_1.5']), axis=1)
    df['integral_error_1.5'] = df.apply(lambda row: _format_error(row['integral_error_1.5']), axis=1)

    _f = io.StringIO()
    df[['pdbID-pdbChain', 'receptor_aa_1-index', 'receptor_aa_value', 'receptor_atom_id', 'receptor_atom_value', 'ligand_id', 'ligand_atom_value', 'euclidean_distance', 'overlap_vdw_radii', 'integral_error_vdw_radii', 'overlap_1.5', 'integral_error_1.5', 'full_receptor_sequence']].to_csv(_f, sep='\t', header=False, index=False)
    _f.seek(0)
    f.write(_f.read().encode('utf8'))
    f.close()


def minf1f2(x, mu2, sd1, sd2):
    min_value = min(e**(-(x**2)/(2*sd1**2))/sd1, e**(-(x-mu2)**2/(2*sd2**2))/sd2)
    return NORM_PDF_FACTOR * min_value


def overlap_radii(dist, r1, r2):
    return quad(
        func=minf1f2,
        a=-np.inf,
        b=np.inf,
        args=(
            dist,
            r1,
            r2
        ),
        epsabs=tol,
        epsrel=tol
    )


def create_distance_file(pdb_id, pdb_chains, receptor_filepaths, ligand_ids, ligand_filepaths, distance_filepath, include_backbone=False, distance_cutoff=20, compressed=True, calculate_overlap=True):

    _rows = []
    previous_pdb_chain = None
    for pdb_chain, receptor_filepath, ligand_id, ligand_filepath in zip(pdb_chains, receptor_filepaths, ligand_ids, ligand_filepaths):

        atoms = defaultdict(list)
        positions = {}

        _last_residue_number = None  # TODO: REMOVE if IGNORE_INSERTION_CODE==False

        receptor_atoms = parsePDB(receptor_filepath, altloc='A1')  # TODO: Find breaking case if altloc is not specified
        receptor_selection = receptor_atoms.select('stdaa and not hydrogen')

        for residue in HierView(receptor_selection).iterResidues():
            residue_number = residue.getResnum()
            residue_name = residue.getResname()

            # positions[residue_number] = AAMAP[residue_name]  # TODO: ADD if IGNORE_INSERTION_CODE==True

            _flag = False  # TODO: REMOVE once the first-4 assumption is removed
            if not include_backbone:
                if _last_residue_number != residue_number:  # TODO: REMOVE if IGNORE_INSERTION_CODE==False
                    positions[residue_number] = AAMAP[residue_name]  # TODO: REMOVE if IGNORE_INSERTION_CODE==False
                    if set(residue.getNames()[:4]) != set(['N', 'CA', 'C', 'O']):  # TODO: REMOVE once the first-4 assumption is removed
                        logger.error(f'Backbone ERROR in {receptor_filepath}-{residue_number}-{residue_name}')
                        _flag = True
                    else:
                        residue = residue.select('not backbone')

            _last_residue_number = residue_number  # TODO: REMOVE if IGNORE_INSERTION_CODE==False
            if residue is not None:
                if _flag:  # TODO: REMOVE once the first-4 assumption is removed
                    atoms[residue_number].extend(
                        [tuple(atom.getCoords()) + (atom.getElement(),) for atom in residue.iterAtoms()][4:]
                    )
                else:
                    atoms[residue_number].extend(
                        [tuple(atom.getCoords()) + (atom.getElement(),) for atom in residue.iterAtoms()]
                    )

        # Sort all atom entries according to distance
        for k, v in atoms.items():
            atoms[k] = sorted(v)

        start_index = min(positions.keys())

        ligand_atoms = parsePDB(ligand_filepath)
        ligand_selection = ligand_atoms.select('hetatm and not hydrogen')

        with tqdm(total=ligand_selection.numAtoms()) as pbar:

            for ligand_atom in ligand_selection:
                ligand_element = ligand_atom.getElement()
                # logger.info(f'Processing ligand {ligand_element}')

                # Special processing for nucleic acids
                if ligand_id == 'NUC':
                    ligand_atom_name = ligand_atom.getName()

                    # Get other atom names in the same chain and residue as *this* atom
                    other_atom_names = HierView(
                        ligand_atom.getAtomGroup()
                    ).getResidue(
                        chid=ligand_atom.getChid(),
                        resnum=ligand_atom.getResnum()
                    ).getNames()

                    # TODO: There has to be a more direct way to do this!
                    ligand_sub_type = 'RNA' if "O2'" in other_atom_names else 'DNA'
                    if "'" in ligand_atom_name or 'P' in ligand_atom_name:
                        ligand_sub_type += 'B'
                else:
                    ligand_sub_type = ''

                for i, ts in atoms.items():
                    n = len(ts)
                    for j, t in enumerate(ts, start=1):
                        _x, _y, _z, _e = t
                        dist = np.linalg.norm(np.array([_x, _y, _z]) - ligand_atom.getCoords())
                        if dist <= distance_cutoff:
                            _row = {
                                'pdbID-pdbChain': f'{pdb_id}{pdb_chain}',
                                'receptor_aa_1-index': i - start_index + 1,
                                'receptor_aa_value': f'{positions[i]}',
                                'receptor_atom_id': f'{j}/{n}',
                                'receptor_atom_value': f'{_e}',
                                'ligand_id': f'{ligand_id}{ligand_sub_type}',
                                'ligand_atom_value': f'{ligand_element}',
                                'euclidean_distance': dist,
                                'full_receptor_sequence': ''
                            }

                            if previous_pdb_chain != pdb_chain:
                               _row['full_receptor_sequence'] = ''.join([positions[aa_index] if aa_index in positions else 'X' for aa_index in range(start_index, max(positions.keys())+1)])
                               previous_pdb_chain = pdb_chain

                            _rows.append(_row)

                pbar.update(1)
    pbar.close()

    df = pd.DataFrame(_rows)

    if calculate_overlap:
        tqdm.pandas()
        logger.info('adding receptor atom vdw radius')
        df['receptor_atom_radius'] = df.progress_apply(lambda row: vdw_radius(row['receptor_atom_value']), axis=1)
        logger.info('adding ligand atom vdw radius')
        df['ligand_atom_radius'] = df.progress_apply(lambda row: vdw_radius(row['ligand_atom_value']), axis=1)

        logger.info('calculating vdw overlap areas')
        df[['overlap_vdw_radii', 'integral_error_vdw_radii']] = df.progress_apply(lambda row: pd.Series(overlap_radii(row['euclidean_distance'], row['receptor_atom_radius'], row['ligand_atom_radius'])), axis=1)
        logger.info('calculating standard overlap areas')
        df[['overlap_1.5', 'integral_error_1.5']] = df.progress_apply(lambda row: pd.Series(overlap_radii(row['euclidean_distance'], 1.5, 1.5)), axis=1)

        if EMULATE_CACHING_BUG:
            overlaps = {}
            logger.info('"caching" calculated overlap areas to emulate legacy behavior')
            for idx in df.index:
                row = df.iloc[idx]
                sd1, sd2 = sorted([vdw_radius(row['receptor_atom_value']), vdw_radius(row['ligand_atom_value'])])
                overlaps[(row['euclidean_distance'], (sd1, sd2))] = row['overlap_vdw_radii'], row['integral_error_vdw_radii']
                overlaps[(row['euclidean_distance'], (1.5, 1.5))] = row['overlap_1.5'], row['integral_error_1.5']

            logger.info('modifying calculated overlap areas to emulate legacy behavior')
            df[['overlap_vdw_radii', 'integral_error_vdw_radii']] = df.progress_apply(lambda row: pd.Series(overlaps[(row['euclidean_distance'], tuple(sorted([row['receptor_atom_radius'], row['ligand_atom_radius']])))]), axis=1)
            df[['overlap_1.5', 'integral_error_1.5']] = df.progress_apply(lambda row: pd.Series(overlaps[(row['euclidean_distance'], (1.5, 1.5))]), axis=1)

    distance_df_to_csv(df, pdb_id, distance_filepath, compressed=compressed)
    return df


def create_fasta(df, filepath, compressed=True, distance='maxstd', distance_cutoff=20):

    distances = {}  # pdbID-pdbChain => <sequence>, [(<residue_index>, <ligand_id>, <distance>), ..]

    for pdb_chain, df1 in df.groupby(['pdbID-pdbChain']):

        receptor_sequence = df1['full_receptor_sequence'].dropna().iloc[0]

        if distance.endswith('std'):
            df1 = df1[df1['integral_error_1.5'] < df1['overlap_1.5']]
        elif distance.endswith('vdw'):
            df1 = df1[df1['integral_error_vdw_radii'] < df1['overlap_vdw_radii']]

        ligand_distances = []
        for residue_index, df2 in df1.groupby(['receptor_aa_1-index']):

            unique_receptor_positions = df2['receptor_atom_id'].unique()
            n_unique_receptor_positions = len(unique_receptor_positions)

            for ligand_id, df3 in df2.groupby(['ligand_id']):

                if distance=='maxstd':
                    d = df3.groupby('receptor_atom_id').sum()['overlap_1.5'].max()
                elif distance=='maxvdw':
                    d = df3.groupby('receptor_atom_id').sum()['overlap_vdw_radii'].max()
                elif distance=='meanstd':
                    d = df3['overlap_1.5'].sum() / n_unique_receptor_positions
                elif distance=='meanvdw':
                    d = df3['overlap_vdw_radii'].sum() / n_unique_receptor_positions
                elif distance=='sumstd':
                    d = df3['overlap_1.5'].sum()
                elif distance=='sumvdw':
                    d = df3['overlap_vdw_radii'].sum()
                elif distance=='mindist':
                    d = df3['euclidean_distance'].min()
                elif distance=='meandist':
                    d = df3.groupby('receptor_atom_id')['euclidean_distance'].min()
                    d = d.reindex(unique_receptor_positions, fill_value=distance_cutoff).mean()
                elif distance=='fracin4':
                    n_aa_atoms = int(df3.iloc[0]['receptor_atom_id'].split('/')[-1])
                    n_aa_atoms_near = (df3.groupby('receptor_atom_id')['euclidean_distance'].min() < 4).sum()
                    d = n_aa_atoms_near/n_aa_atoms

                ligand_distances.append((residue_index, ligand_id, d))

        distances[pdb_chain] = receptor_sequence, sorted(ligand_distances)

    with open(filepath, 'w') as f:
        for chain, (seq, ligand_distances) in distances.items():
            f.write(f'>{chain} bindingSiteRes=')
            f.write(','.join([f'{pos}-{lig}-{d:.5f}' for pos, lig, d in ligand_distances]))
            f.write(f';\n{seq}\n\n')
