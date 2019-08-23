import os
import logging
import requests
import json
from subprocess import run
from tempfile import NamedTemporaryFile

logger = logging.getLogger(__name__)


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

    response = requests.post(
        'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan',
        headers={'Content-type': 'application/json', 'Accept': 'application/json'},
        data=json.dumps({
            'hmmdb': 'pfam',
            'cut_ga': True,
            'seq': sequence
        }),
    )

    d = response.json()
    logger.info('Hmmr Web API response obtained')

    hits = d['results']['hits']
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
