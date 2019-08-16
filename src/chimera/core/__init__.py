import logging
import numpy as np
import pandas as pd

from chimera import binding_frequencies_interacdome, binding_frequencies_dsprint
from chimera.core.hmmr import find_hmmr_domains_web

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


def query(sequence, algorithm='dSPRINT'):
    """
    Find out binding frequency data suitable for display on the site
    :param sequence: String of Protein sequence of length L
    :param algorithm: One of 'dSPRINT' or 'InteracDome' (case-sensitive)
    :return: A 2-tuple of values
        0: DataFrame containing Hmmer matches
        1: DataFrame containing Ligand-binding frequency data
    """
    hits = find_hmmr_domains_web(sequence)

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

    domains = pd.DataFrame(results)

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

    binding_frequencies = binding_frequencies_interacdome if algorithm == 'InteracDome' else binding_frequencies_dsprint
    df = pd.merge(matches, binding_frequencies, left_on=['pfam_domain', 'match_i'], right_on=['pfam_id', 'match_state'])
    logger.info('Merged domain results / binding frequency dataFrames')

    return domains, df