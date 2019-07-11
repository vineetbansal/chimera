from collections import defaultdict
import numpy as np
from math import e, pi, sqrt
from scipy.integrate import quad
import gzip
import logging
from tqdm import tqdm

from prody import parsePDB
from prody.atomic.hierview import HierView
from prody.atomic import AAMAP

from chimera.utils import vdw_radius

tol = np.finfo(float).eps ** 0.25
logger = logging.getLogger(__name__)

# TODO: Whether to ignore residue insertion code (col 27 in pdb) when identifying unique residues
# Should ideally be False, but set to True for backwards compatibility
IGNORE_INSERTION_CODE = True

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


def _format_integral(x):
    s = f'{x:.7g}'
    if s=='0.0001':
        return '1e-04'
    else:
        return s


def _format_error(x):
    s = f'{x:.2g}'
    if s=='0.0001':
        return '1e-04'
    else:
        return s


def minf1f2(x, mu2, sd1, sd2):
    min_value = min(e**(-(x**2)/(2*sd1**2))/sd1, e**(-(x-mu2)**2/(2*sd2**2))/sd2)
    return min_value/sqrt(2*pi)


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

    _open = open
    if compressed:
        _open = gzip.open
        if not distance_filepath.endswith('.gz'):
            distance_filepath += '.gz'

    f = _open(distance_filepath, 'wb')
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

                            r1 = e1 = r2 = e2 = 'N/A'
                            if calculate_overlap:
                                _r1, _r2 = vdw_radius(_e), vdw_radius(ligand_element)
                                r1, e1 = overlap_radii(dist, _r1, _r2)
                                r2, e2 = overlap_radii(dist, 1.5, 1.5)

                            f.write(('\t'.join([
                                f'{pdb_id}{pdb_chain}',
                                f'{i - start_index + 1}',
                                f'{positions[i]}',
                                f'{j}/{n}',
                                f'{_e}',
                                f'{ligand_id}{ligand_sub_type}',
                                f'{ligand_element}',
                                f'{dist}',
                                _format_integral(r1),
                                _format_error(e1),
                                _format_integral(r2),
                                _format_error(e2)
                            ])).encode('utf8'))

                            if previous_pdb_chain != pdb_chain:
                                f.write(
                                    ('\t' + ''.join([positions[aa_index] if aa_index in positions else 'X' for aa_index in range(start_index, max(positions.keys())+1)]) +
                                     '\n').encode('utf8'))
                                previous_pdb_chain = pdb_chain
                            else:
                                f.write(b'\t\n')
                pbar.update(1)

    pbar.close()
    f.close()
