import os
import logging
from subprocess import run
from tempfile import NamedTemporaryFile
from Bio import SearchIO

from . import DomainFinder
from chimera import config

logger = logging.getLogger(__name__)


class HmmerDomainFinder(DomainFinder):
    def find_domains(self, sequences):

        f_in = NamedTemporaryFile(delete=False)
        f_alignments = NamedTemporaryFile(delete=False)
        f_out = NamedTemporaryFile(delete=False)

        for sequence in sequences:
            f_in.write(f'>{sequence.name}\n{sequence.seq}\n\n'.encode('utf8'))
        f_in.close()
        logger.info(f'Wrote sequences to temporary file = {f_in.name}')

        pfam_hmm_path = os.path.join(config.dirs.data.pfam, 'Pfam-A.hmm')
        hmmscan_bin = os.path.join(config.dirs.bin.hmmr, 'hmmscan')

        try:
            p = run([
                hmmscan_bin,
                '--cut_ga',
                '-o',
                f_alignments.name,
                '--domtblout',
                f_out.name,
                pfam_hmm_path,
                f_in.name
            ], capture_output=True)

            if p.returncode != 0:
                raise RuntimeError('Hmmscan execution failed')

            logger.info(f'Wrote hmmscan stdout to temporary file = {f_alignments.name}')

            """
            From http://biopython.org/DIST/docs/api/Bio.SearchIO-module.html:

            SearchIO parses a search output file's contents into a hierarchy of four nested objects:
                QueryResult represents a search query. This is the main object returned by the input functions and it
                    contains all other objects.
                Hit represents a database hit
                HSP represents high-scoring alignment region(s) in the hit
                HSPFragment represents a contiguous alignment within the HSP

            The 'fragments' attribute of a QueryResult gives us all HSPFragment objects as a list.
            """
            qresults = SearchIO.parse(f_out.name, format='hmmscan3-domtab')
            alignment_qresults = SearchIO.parse(f_alignments.name, format='hmmer3-text')

            # Create a dictionary of sequence strings so that we can get alignment information of any 'hit' in
            # qresult (corresponding to Hmmer --domtab result) by consulting the corresponding entry in qresult
            # alignments (corresponding to Hmmer stdout)
            # TODO: There's probably a more direct way to do this using the Bio.SearchIO library

            # A dict mapping (query_id, hit_id, hit_start, hit_end, query_start, query_end) => <query_seq_string>
            # Note that <query_seq_string> may contain lowercase letters (for insertions) and hyphens for deletions.
            alignments = {}
            for alignment_qresult in alignment_qresults:
                for f in alignment_qresult.fragments:
                    key = alignment_qresult.id, f.hit_id, f.hit_start, f.hit_end, f.query_start, f.query_end
                    if key in alignments:
                        raise RuntimeError('Unexpected - a key with identical hit info already exists!')
                    alignments[key] = str(f.query.seq)

            results = []
            for qresult in qresults:
                for hit in qresult.hits:
                    for hsp in hit.hsps:
                        result = dict()
                        result['query_id'] = qresult.id
                        result['alihmmacc'] = hit.accession
                        result['alihmmname'] = hit.id
                        result['alisqfrom'] = hsp.query_start + 1  # 0-indexed inclusive -> 1-indexed inclusive
                        result['alisqto'] = hsp.query_end          # 0-indexed exclusive -> 1-indexed inclusive
                        result['alihmmfrom'] = hsp.hit_start + 1   # 0-indexed inclusive -> 1-indexed inclusive
                        result['alihmmto'] = hsp.hit_end           # 0-indexed exclusive -> 1-indexed inclusive
                        result['aliM'] = hit.seq_len
                        result['bitscore'] = hsp.bitscore
                        result['ievalue'] = hsp.evalue

                        # Get alignment info from the dictionary we constructed from qresult_alignments
                        key = qresult.id, hit.id, hsp.hit_start, hsp.hit_end, hsp.query_start, hsp.query_end
                        if key not in alignments:
                            raise RuntimeError(f'Alignment information for {key} not found!')
                        result['aliaseq'] = alignments[key]

                        # TODO: The old Interacdome website does the following filtering of results for Hmmer
                        if result['bitscore'] > 0:
                            results.append(result)

        finally:
            os.unlink(f_in.name)
            os.unlink(f_alignments.name)
            os.unlink(f_out.name)

        return results
