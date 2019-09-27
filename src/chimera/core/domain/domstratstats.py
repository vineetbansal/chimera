import os
import logging
from subprocess import run
from tempfile import NamedTemporaryFile
from collections import defaultdict
from Bio import SearchIO

from chimera.core.domain import DomainFinder
from chimera import config

logger = logging.getLogger(__name__)


class DomStratStatsDomainFinder(DomainFinder):
    def find_domains(self, sequences):
        def read_filtered_domtab(filename):
            """
            Read a 'filtered' Hmmer domtab file (i.e. a Hmmer domtab files with arbitrary lines removed) as a dictionary.
            This is the kind of file that DomStratStats produces as output.

            Unfortunately, Bio.SearchIO is incapable of parsing such a file since it relies too heavily on the lines in the
            domtab file following a certain order (lines sorted by query name, then by target name, then by domain index,
            with no gaps in between).

            Here, we simply create a dictionary of the following form:

                {<query_id> => { <hit_id> => [<domain_index>, <domain_index>, ..], .. } .. }

            where <query_id> is the query_name (e.g. 'ctcf'), and <hit_id> is the target name (e.g. 'zf-H2C2_2')
            and <domain_index> is the hit # of this particular domain for this query sequence (e.g. 3) that dPUC2 deems
            important enough to include in it's results.

            :param filename: output file produced by dPUC2 (by 1dpuc2.pl script)
            :return: A dictionary that follows the following structure:
                {<query_id> => { <hit_id> => [<domain_index>, <domain_index>, ..], .. } .. }
            """

            # https://stackoverflow.com/questions/5029934/python-defaultdict-of-defaultdict
            query_dict = defaultdict(lambda: defaultdict(list))

            if not os.path.exists(filename):
                filename = filename + '.gz'
                if not os.path.exists(filename):
                    raise RuntimeError(f'Input file {filename} not found!')

            from smart_open import open

            with open(filename, 'r') as f:
                for line in f.readlines():
                    line = line.strip()
                    if not line.startswith('#'):
                        cols = [x for x in line.strip().split(' ') if x]
                        hit_id = cols[0]
                        query_id = cols[3]
                        domain_index = int(cols[9])

                        query_dict[query_id][hit_id].append(domain_index)

            return query_dict

        # Since DomStratStats insists on appending '.gz' to filenames and refuses to work
        # if the corresponding filename (without .gz) already exists, we execute the following statements just to get
        # some temporary filenames, by creating them and then removing them immediately.
        # TODO: delete=True can be switched to delete=False once DomStratStats supports a '--noGzip' option like dPuc2
        f_in = NamedTemporaryFile(delete=False)
        f_intermediate1 = NamedTemporaryFile(delete=True)
        f_intermediate1.close()
        f_intermediate2 = NamedTemporaryFile(delete=True)
        f_intermediate2.close()
        f_out = NamedTemporaryFile(delete=True)
        f_out.close()

        f_alignments = NamedTemporaryFile(delete=False)

        for sequence in sequences:
            f_in.write(f'>{sequence.name}\n{sequence.seq}\n\n'.encode('utf8'))
        f_in.close()
        logger.info(f'Wrote sequences to temporary file = {f_in.name}')

        pfam_hmm_path = os.path.join(config.dirs.data.pfam, 'Pfam-A.hmm')
        pfam_hmm_dat_path = os.path.join(config.dirs.data.pfam, 'Pfam-A.hmm.dat')
        hmmscan_bin = os.path.join(config.dirs.bin.hmmr, 'hmmscan')

        try:
            p1 = run([
                'perl',
                '0runHmmscan.pl',
                hmmscan_bin,
                pfam_hmm_path,
                f_in.name,
                f_intermediate1.name,
                f_alignments.name
            ], cwd=config.dirs.bin.domstratstats, capture_output=True)

            if p1.returncode != 0:
                raise RuntimeError('Perl execution failed')

            logger.info(f'Wrote DomStratStats hmmscan results to temporary file = {f_intermediate1.name}')
            logger.info(f'Wrote hmmscan alignments to temporary file = {f_alignments.name}')

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

            p2 = run([
                'perl',
                '1noOvs.pl',
                f_intermediate1.name,
                f_intermediate2.name,
                pfam_hmm_dat_path
            ], cwd=config.dirs.bin.domstratstats, capture_output=True)

            if p2.returncode != 0:
                raise RuntimeError('Perl execution failed')

            logger.info(f'Wrote DomStratStats results to output file = {f_intermediate2.name}')

            p3 = run([
                'perl',
                '2domStratStats.pl',
                f_in.name,
                f_intermediate2.name,
                f_out.name
            ], cwd=config.dirs.bin.domstratstats, capture_output=True)

            if p3.returncode != 0:
                raise RuntimeError('Perl execution failed')

            logger.info(f'Wrote DomStratStats results to output file = {f_out.name}')

            logger.info(f'Obtaining list of dPUC2 hits from output file.')
            final_results = read_filtered_domtab(f_out.name)

            results = []

            qresults = SearchIO.parse(f_intermediate1.name, format='hmmscan3-domtab')
            for qresult in qresults:
                for hit in qresult.hits:
                    for hsp in hit.hsps:

                        # Is this HSP included in dPUC2 results? If so, add it to our results
                        if hsp.domain_index in final_results[qresult.id][hit.id]:
                            result = dict()
                            result['query_id'] = qresult.id
                            result['alihmmacc'] = hit.accession
                            result['alihmmname'] = hit.id
                            result['alisqfrom'] = hsp.query_start + 1  # 0-indexed inclusive -> 1-indexed inclusive
                            result['alisqto'] = hsp.query_end  # 0-indexed exclusive -> 1-indexed inclusive
                            result['alihmmfrom'] = hsp.hit_start + 1  # 0-indexed inclusive -> 1-indexed inclusive
                            result['alihmmto'] = hsp.hit_end  # 0-indexed exclusive -> 1-indexed inclusive
                            result['aliM'] = hit.seq_len
                            result['bitscore'] = hsp.bitscore
                            result['ievalue'] = hsp.evalue

                            # Get alignment info from the dictionary we constructed from qresult_alignments
                            key = qresult.id, hit.id, hsp.hit_start, hsp.hit_end, hsp.query_start, hsp.query_end
                            if key not in alignments:
                                raise RuntimeError(f'Alignment information for {key} not found!')
                            result['aliaseq'] = alignments[key]

                            results.append(result)

        finally:
            os.unlink(f_in.name)
            os.unlink(f_alignments.name)

            # TODO: This logic can be simplified once DomStratStats supports a '--noGzip' option
            for f in f_intermediate1, f_intermediate2, f_out:
                fname = f.name
                if os.path.exists(fname):
                    os.unlink(fname)
                if os.path.exists(fname + '.gz'):
                    os.unlink(fname + '.gz')

        return results
