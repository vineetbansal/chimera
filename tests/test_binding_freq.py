from unittest import TestCase
from importlib.resources import read_text

from chimera.utils import parse_fasta
from chimera.core import seq_to_matchstates

ctcf = read_text('chimera.data.sample', 'ctcf.fa')
ctcf = str(parse_fasta(ctcf)[0].seq)


class BindingSeqTestCase(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testMatchStatesNormal(self):
        # These are actual domain hits from Hmmr for CTCF, and represent aliseq/target_start/target_end respectively
        domain_seq, start, end = 'FQCELCSYTCPRRSNLDRHMKSH', 266, 288  # domain_seq length = 23

        # start/end indices are both included in the domain-hit, so the following holds
        self.assertEqual(len(domain_seq), end-start+1)

        # start/end are 1-indexed - get the corresponding fragment from the original (CTCF) sequence
        ctcf_fragment = ctcf[start-1: end]

        # In this case, this is an exact match
        self.assertEqual(domain_seq, ctcf_fragment)

        match_states, seq_indices = seq_to_matchstates(domain_seq, start, end)

        # match_states are contiguous, 1..23
        self.assertEqual(match_states, tuple(range(1, 24)))

    def testMatchStatesInsertion(self):
        # lowercase letters represent insertions wrt the domain, in this case 'k' at 1-indexed position 23
        domain_seq, start, end = 'HKCPDCDMAFVTSGELVRHRRYkH', 322, 345  # domain_seq length = 24

        # start/end indices are both included in the domain-hit, so the following still holds
        self.assertEqual(len(domain_seq), end - start + 1)

        # start/end are 1-indexed - get the corresponding fragment from the original (CTCF) sequence
        ctcf_fragment = ctcf[start - 1: end]

        # In this case, this is NOT an exact match
        self.assertNotEqual(domain_seq, ctcf_fragment)
        # But IS an exact match if we convert the insertions found by Hmmer to uppercase
        self.assertEqual(domain_seq.upper(), ctcf_fragment)

        match_states, seq_indices = seq_to_matchstates(domain_seq, start, end)

        # match_states are contiguous, 1..23 (insertions wrt the domain are ignored)
        self.assertEqual(match_states, tuple(range(1, 24)))

    def testMatchStatesDeletion(self):
        # Hyphens represent deletions wrt the domain, in this case the '-' at 1-indexed position 22
        domain_seq, start, end = 'FQCSLCSYASRDTYKLKRHMR-THSG', 379, 403  # domain_seq length = 26

        # Because of the deletion, the ctcf fragment is 1 short in length wrt the domain sequence
        self.assertEqual(len(domain_seq) - 1, end - start + 1)

        # start/end are 1-indexed, get the corresponding fragment from the original (CTCF) sequence
        ctcf_fragment = ctcf[start - 1: end]

        # In this case, this is NOT an exact match
        self.assertNotEqual(domain_seq, ctcf_fragment)
        # But IS an exact match if we strip out the '-' from the domain
        self.assertEqual(domain_seq.replace('-', ''), ctcf_fragment)

        match_states, seq_indices = seq_to_matchstates(domain_seq, start, end)

        # match_states are NOT contiguous, 1..21, 23..26 (i.e. missing 22, the 1-indexed position of the deletion)
        self.assertEqual(match_states, tuple(range(1, 22)) + tuple(range(23, 27)))
