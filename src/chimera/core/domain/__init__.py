import pandas as pd


class DomainFinder:
    """
    A class that finds domains, given one of more sequences.
    """

    def find_domains(self, sequences):
        """
        Find domains for one or more sequences
        :param sequence: A list of Bio.SeqRecord.SeqRecord objects
        :return: A list of dicts with keys
            name: str (name of domain), e.g. "zf-H2C2_2"
            ndom: int (no. of instances of this domain found)
            domains: A list of dicts with keys:
                query_id - str, query id, e.g. 'ctcf'
                alihmmacc - accession (str), e.g. 'PF13465.6'
                alihmmname - str, e.g. 'zf-H2C2_2'
                alisqfrom - int. e.g. 262
                alisqto - int, e.g. 274
                alihmmfrom - int, e.g. 11
                alihmmto - int, e.g. 23
                aliM - int, e.g. 26
                bitscore - float, e.g. 5.896
                ievalue - str, e.g. '19'
                aliaseq - str, e.g. 'VKKTFQCELCSYT'
        """
        raise NotImplementedError('Subclasses must implement this.')

    def domain_table(self, sequences, full_domains=False):

        hits = self.find_domains(sequences)
        results = []  # A list-of-dicts that we'll convert to a DataFrame
        for hit in hits:
            if full_domains and not (hit['alihmmfrom'] == 1 and hit['alihmmto'] == hit['aliM']):
                continue
            pfam_name = hit['alihmmacc'][:7] + '_' + hit['alihmmname']  # TODO: Why this strange clipping of names?
            results.append({
                'query_id': hit['query_id'],
                'pfam_domain': pfam_name,
                'target_start': int(hit['alisqfrom']),
                'target_end': int(hit['alisqto']),
                'hmm_start': int(hit['alihmmfrom']),
                'hmm_end': int(hit['alihmmto']),
                'domain_length': int(hit['aliM']),
                'bit_score': float(hit['bitscore']),
                'e_value': float(hit['ievalue']),
                'aliseq': hit['aliaseq'],
            })

        return pd.DataFrame(results)


from . hmmer import HmmerDomainFinder
from . hmmerweb import HmmerWebDomainFinder
from . dpuc2 import Dpuc2DomainFinder
from . domstratstats import DomStratStatsDomainFinder
