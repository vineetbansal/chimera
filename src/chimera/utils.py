import requests
import json
import logging
import importlib.resources
import pandas as pd
import os
from subprocess import run
from io import StringIO
from tempfile import NamedTemporaryFile
from Bio import SeqIO

import chimera.data

logger = logging.getLogger(__name__)

vdw_df = None
with importlib.resources.path(chimera.data, 'vdw.txt') as path:
    vdw_df = pd.read_csv(
        path,
        sep='\t',
        header=None,
        comment='#',
        names=(
            'atomic_number',
            'element_symbol',
            'atomic_radius',
            'ionic_radius',
            'covalent_radius',
            'vdw_radius',
            'crystal_radius'
        ),
        index_col=1
    )
    vdw_df.index = vdw_df.index.str.lower()
    vdw_df = vdw_df.replace({pd.np.nan: None})


def vdw_radius(element):
    """
    Return VDW radius of element
    :param element: Element symbol (case-sensitive)
    :return: Float indicating the VDW radius
    """
    s = vdw_df.loc[element.lower()]
    return s.vdw_radius or s.ionic_radius or 1.5


_ligand_to_groups = {}


def ligand_to_groups(ligand_id):

    global _ligand_to_groups
    if not _ligand_to_groups:
        with importlib.resources.path(chimera.data, 'ligand_groups.txt') as path:
            df = pd.read_csv(
                path,
                sep='\t',
                header=None,
                comment='#',
                names=('group_name', 'ligand_id'),
                na_filter=False  # important since 'NA' is a valid ligand identifier!
            )

            _ligand_to_groups = {}
            for l_id, _df in df.groupby('ligand_id'):
                group_names = set(_df.group_name.values)

                # a single molecule must belong exclusively to the nucleic acid, ion, or small molecule groups
                if 'NUCACID_' in group_names:
                    [group_names.discard(g) for g in ['ION_', 'METABOLITE_', 'DRUGLIKE_', 'SM_']]
                elif 'ION_' in group_names:
                    [group_names.discard(g) for g in ['METABOLITE_', 'DRUGLIKE_', 'SM_']]

                _ligand_to_groups[l_id] = sorted(group_names)

    group_names = ['ALL_', ligand_id] + _ligand_to_groups[ligand_id]
    if not ('NUCACID_' in group_names or 'ION_' in group_names or 'III' in group_names):
        group_names.append('SM_')

    # Last-minute Ligand Group renamings
    group_renamings = {
        'NUCDNA': 'DNABASE_',
        'NUCDNAB': 'DNABACKBONE_',
        'NUCRNA': 'RNABASE_',
        'NUCRNAB': 'RNABACKBONE_',
        'III': 'PEPTIDE_'
    }

    group_names = [group_renamings.get(g, g) for g in group_names]

    return group_names


def find_hmmr_domains_web(sequence):
    """
    Search a sequence against Hmmr Pfam Profile Database by making a Hmmer Web API call
    :param sequence: protein chain sequence, a string
    :return: A list of dicts with keys
        name (name of domain)
        ndom (no. of instances of this domain found)
        domains: A list of dicts with keys:
            alihmmacc
            alihmmname
            alisqfrom
            alisqto
            alihmmfrom
            alihmmto
            aliM
            bitscore
            is_reported
            ievalue
            aliaseq
    """

    r = requests.post(
        'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan',
        headers={'Content-type': 'application/json', 'Accept': 'application/json'},
        data=json.dumps({
            'hmmdb': 'pfam',
            'cut_ga': True,
            'seq': sequence
        }),
    ).json()
    logger.info('Hmmr Web API response obtained')

    hits = r['results']['hits']
    return hits


def find_hmmr_domains_local(sequence):
    """
    Search a sequence against Hmmr Pfam Profile Database by making a local call to Hmmr
    :param sequence: protein chain sequence, a string
    :return: A list of dicts with keys
        name (name of domain)
        ndom (no. of instances of this domain found)
        domains: A list of dicts with keys:
            alihmmacc
            alihmmname
            alisqfrom
            alisqto
            alihmmfrom
            alihmmto
            aliM
            bitscore
            is_reported
            ievalue
            aliaseq
    """

    # TODO: Abstract out
    pfam_path = '/home/vineetb/git_checkouts/run-hmmer/Pfam-A.hmm'
    f = NamedTemporaryFile(delete=False)
    f.write(sequence.encode('utf8'))
    f.close()

    p = run(['/opt/hmmr/bin/hmmscan', pfam_path, f.name], capture_output=True)
    s = p.stdout.decode('utf8')
    print(p.returncode)

    os.unlink(f.name)


def parse_fasta(s, default_seq_id='seq0'):
    sequences = list(SeqIO.parse(StringIO(s), 'fasta'))
    if not sequences:
        sequences = list(SeqIO.parse(StringIO(f'>{default_seq_id}\n' + s), 'fasta'))
    return sequences
