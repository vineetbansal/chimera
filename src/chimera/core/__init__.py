import logging
import numpy as np
import pandas as pd

from chimera import binding_frequencies_interacdome, binding_frequencies_dsprint
from chimera.core.domain import HmmerDomainFinder, HmmerWebDomainFinder, Dpuc2DomainFinder

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


def query(sequences, algorithm='dsprint', domain_algorithm='hmmer'):
    """
    Find out binding frequency data suitable for display on the site
    :param sequences: A list of Bio.SeqRecord.SeqRecord objects
    :param algorithm: Ligand-binding algorithm (case-sensitive). One of:
        dsprint
        interacdome
    :param domain_algorithm: Domain-finding algorithm (case-sensitive). One of:
        hmmer
        hmmerweb
        dpuc2
    :return: A 2-tuple of values
        0: DataFrame containing Hmmer matches
        1: DataFrame containing Ligand-binding frequency data
    """

    domain_finder = {
        'hmmer': HmmerDomainFinder,
        'hmmerweb': HmmerWebDomainFinder,
        'dpuc2': Dpuc2DomainFinder
    }[domain_algorithm]()

    domains = domain_finder.domain_table(sequences)

    domains[['match_states', 'seq_indices']] = domains.apply(
        lambda row: pd.Series(seq_to_matchstates(row.aliseq, row.target_start, row.target_end)),
        axis=1
    )
    logger.info('Added match state and sequence index information to results')

    match_rows = []
    for (query_id, pfam_domain), _df in domains.groupby(['query_id', 'pfam_domain']):
        for _, row in _df.iterrows():
            for match_i, seq_i in zip(row.match_states, row.seq_indices):
                match_rows.append({'query_id': query_id, 'pfam_domain': pfam_domain, 'match_i': match_i, 'seq_i': seq_i})
    matches = pd.DataFrame(match_rows)
    logger.info('Created an unpivoted match/sequence table of domain results')

    if algorithm == 'interacdome':
        binding_frequencies = binding_frequencies_interacdome
    elif algorithm == 'dsprint':
        binding_frequencies = binding_frequencies_dsprint
    else:
        raise RuntimeError('Unsupported ligand frequency algorithm')

    df = pd.merge(matches, binding_frequencies, left_on=['pfam_domain', 'match_i'], right_on=['pfam_id', 'match_state'])
    logger.info('Merged domain results / binding frequency dataFrames')

    return domains, df
