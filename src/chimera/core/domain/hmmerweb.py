import logging
import requests
import json

from . import DomainFinder

logger = logging.getLogger(__name__)


class HmmerWebDomainFinder(DomainFinder):

    def find_domains(self, sequences):
        results = []
        for sequence in sequences:

            # TODO: Does HmmerWeb support a single POST call with multiple sequences?
            logger.info(f'Making Hmmer web call for sequence {sequence.name}')
            response = requests.post(
                'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan',
                headers={'Content-type': 'application/json', 'Accept': 'application/json'},
                data=json.dumps({
                    'hmmdb': 'pfam',
                    'cut_ga': True,
                    'seq': str(sequence.seq)
                }),
            )

            d = response.json()
            logger.info(f'Hmmr Web API response for sequence {sequence.name} obtained')

            hits = d['results']['hits']

            # The Hmmer Web call returns way more results than a locally run 'hmmscan --cut_ga'.
            # This is because it includes records with 'is_included' = 0, which we need to filter out

            # TODO: There's probably a way to enhance the POST call above with extra params, to directly take
            # care of this.

            for hit in hits:

                # 'nincluded' key for each 'hit' tells us how many records in the 'domains' key
                # will have is_included = 1, so we need not enter the inner loop

                if hit['nincluded'] > 0:
                    for record in hit['domains']:
                        if record['is_included']:
                            result = {
                                'query_id': sequence.name,
                                'alihmmacc': record['alihmmacc'],
                                'alihmmname': record['alihmmname'],
                                'alisqfrom': int(record['alisqfrom']),
                                'alisqto': int(record['alisqto']),
                                'alihmmfrom': int(record['alihmmfrom']),
                                'alihmmto': int(record['alihmmto']),
                                'aliM': int(record['aliM']),
                                'bitscore': float(record['bitscore']),
                                'ievalue': float(record['ievalue']),
                                'aliaseq': record['aliaseq'],
                            }

                            # TODO: The old Interacdome website does the following filtering of results for Hmmer
                            # This does not seem applicable while running dPuc2
                            if result['bitscore'] > 0 and result['alihmmfrom'] == 1 and \
                                    result['alihmmto'] == result['aliM']:
                                results.append(result)

        return results
