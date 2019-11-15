from tempfile import NamedTemporaryFile
from celery import Celery

from chimera import config
from chimera.core import query as sequence_query
from chimera.utils import email, parse_fasta

app = Celery('tasks', broker=config.celery.broker, backend=config.celery.backend)
app.control.enable_events()


@app.task
def query(sequences=None, seq_text=None, save_results=False, email_address=None, algorithm='dsprint', domain_algorithm='hmmer', silent=False):
    """
    Query 1 or more sequences and send final results to an email address

    :param sequences: A list of Bio.SeqRecord.SeqRecord objects
    :param seq_text: Text corresponding to one or more sequences. Takes precedence over 'sequences'.
        Useful when running this function 'delayed' to avoid serializing the 'sequences' array.
    :param save_results: Whether to save results in a .csv file
    :param email_address: Recipient's email address where we send results once processing is complete
    :param algorithm: Ligand-binding algorithm (case-sensitive). One of:
        dsprint
        interacdome
    :param domain_algorithm: Domain-finding algorithm (case-sensitive). One of:
        hmmer
        hmmerweb
        dpuc2
    :param silent: A boolean indicating whether we return anything
        Useful when this function is executed in asynchronous mode to avoid serialization issues with returned objects
    :return: If silent is False, a 3-tuple of values
        0: DataFrame containing Hmmer matches
        1: DataFrame containing Ligand-binding frequency data
        2: The file path of saved results if save_results=True, None otherwise
    """
    if seq_text is not None:
        sequences = parse_fasta(seq_text)
    sequences = sequences[:config.web.max_sequences]

    if email_address is not None and not save_results:
        raise RuntimeError('save_results should be True if an email address is specified.')

    domains, df = sequence_query(sequences, algorithm=algorithm, domain_algorithm=domain_algorithm)

    result_filename = None
    if save_results:
        with NamedTemporaryFile(delete=False) as f:
            domains.to_csv(f.name, index=False)
        result_filename = f.name

    if email_address is not None and email_address.strip():
        email(email_address, f'ProtDomain Results for job {query.request.id}', 'ProtDomain Results are attached.', [result_filename])

    if not silent:
        return domains, df, result_filename

