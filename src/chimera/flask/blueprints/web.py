import logging
import numpy as np
import pandas as pd
import requests
import json
from flask import Blueprint, request, render_template

from chimera import config, binding_frequencies_interacdome, binding_frequencies_dsprint
from chimera.data.sample import ctcf
from chimera.utils import find_hmmr_domains_web
from chimera.plots import binding_freq_plot_data_domain, binding_freq_plot_data_sequence

bp = Blueprint('web', __name__)
logger = logging.getLogger(__name__)


def seq_to_matchstates(seq, start, end):
    """
    Determine the 'index' and 'matchstate' information of a given sequence.
    TODO: Code ported directly from R - Could use an intuitive explanation!
    :param seq: A string of single-letter components of a sequence, which can contain lowercase characters
        or '-' characters.
    :return: A 2-tuple of arrays
        0: array indicating indices
        1: array indicating match-states.
    """
    match_states = np.zeros(len(seq)).astype('int')
    match_states_mask = np.array([s_i.upper() == s_i for s_i in seq])
    match_states[match_states_mask] = np.arange(1, len(match_states[match_states_mask])+1)
    seq_indices = np.zeros(len(seq)).astype('int')
    seq_indices_mask = np.array([s_i != '-' for s_i in seq])
    seq_indices[seq_indices_mask] = np.arange(start, end+1)
    assert len(match_states) == len(seq_indices), "Match states and seq indices must have same length"

    # line up the outputs, filter out any 0s (in either of them)
    match_states, seq_indices = zip(*[(m, s) for m, s in zip(match_states, seq_indices) if m != 0 and s != 0])

    return match_states, seq_indices


@bp.route('/interacdome', methods=['GET', 'POST'])
def interacdome():
    from chimera import df_dl

    df_dl = df_dl[
        (df_dl.num_nonidentical_instances >= config.web.min_instances) &
        (df_dl.num_structures >= config.web.min_structures)
    ]

    pfam_ids = pd.unique(df_dl['pfam_id'])

    selected_pfam_id = None
    data = []
    if request.method == 'POST':
        selected_pfam_id = request.form['pfam_id']
        data = binding_freq_plot_data_domain(selected_pfam_id)

    return render_template('interacdome.html', pfam_ids=pfam_ids, selected_pfam_id=selected_pfam_id, data=data)


@bp.route('/faqs')
def faqs():
    return render_template('faqs.html')


def _query(sequence, algorithm='dsprint'):
    """
    Find out binding frequency data suitable for display on the site
    :param sequence: String of Protein sequence of length L
    :param algorithm: One of 'dsprint' or 'interacdome'
    :return: A 2-tuple of values
        0: Iterable of strings, indicating ligand-types (length M)
        1: An M x L ndarray of floats, indicating binding frequencies corresponding to each ligand-type/position
            combination.
    """
    sequence_length = len(sequence)
    hits = find_hmmr_domains_web(sequence)['hits']

    results = []  # A list-of-dicts that we'll convert to a DataFrame
    for hit in hits:
        for d in hit['domains']:
            pfam_name = d['alihmmacc'][:7] + '_' + d['alihmmname']  # TODO: Why this strange clipping of names?
            results.append({
                'pfam_domain': pfam_name,
                'target_start': int(d['alisqfrom']),
                'target_end': int(d['alisqto']),
                'hmm_start': int(d['alihmmfrom']),
                'hmm_end': int(d['alihmmto']),
                'domain_length': int(d['aliM']),
                'bit_score': float(d['bitscore']),
                'reported': bool(d['is_reported']),
                'e_value': float(d['ievalue']),
                'aliseq': d['aliaseq'],
            })

    return pd.DataFrame(results)

    # TODO - if present in df_bp['pfam_id'].unique()
    # TODO - Should be same as one in interacdome_allresults (and thus binding_frequencies_*.csv)
    domains['interacdome'] = False

    domains = domains[(
        domains['bit_score'] > 0) &
        (domains['hmm_start'] == 1) &
        (domains['hmm_end'] == domains['domain_length'])
    ]
    logger.info(f'After filtering, obtained {len(domains)} domain results from Hmmr.')

    domains[['match_states', 'seq_indices']] = domains.apply(
        lambda row: pd.Series(seq_to_matchstates(row.aliseq, row.target_start, row.target_end)),
        axis=1
    )
    logger.info('Added match state and sequence index information to results')

    match_rows = []
    for pfam_domain, _df in domains.groupby('pfam_domain'):
        for _, row in _df.iterrows():
            for match_i, seq_i in zip(row.match_states, row.seq_indices):
                match_rows.append({'pfam_domain': pfam_domain, 'match_i': match_i, 'seq_i': seq_i})
    matches = pd.DataFrame(match_rows)
    logger.info('Created an unpivoted match/sequence table of domain results')

    binding_frequencies = binding_frequencies_interacdome if algorithm == 'interacdome' else binding_frequencies_dsprint
    df = pd.merge(matches, binding_frequencies, left_on=['pfam_domain', 'match_i'], right_on=['pfam_id', 'match_state'])
    logger.info('Merged domain results / binding frequency dataFrames')

    ligand_types = df.ligand_type.unique()

    data = np.zeros((len(ligand_types), sequence_length))
    for i, ligand_type in enumerate(ligand_types):
        bf = df[df.ligand_type == ligand_type].groupby('seq_i')['binding_frequency'].max()
        bf = bf.reindex(pd.RangeIndex(1, sequence_length+1), fill_value=0)
        data[i, :] = bf.values

    # TODO: When displaying, sort by target_start (df = df.sort_values('target_start'))

    return ligand_types, data


@bp.route('/', methods=['GET', 'POST'])
def index():
    data = []
    if request.method == 'POST':
        seq = request.form['seqTextArea']
        algorithm = request.form['algorithmSelect'].lower()
        ligand_types, data = _query(seq, algorithm)
        data = binding_freq_plot_data_sequence(seq, ligand_types, data, algorithm)

    return render_template('index.html', data=data, sample_seq=ctcf)
