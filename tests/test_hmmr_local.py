from unittest import TestCase
from importlib.resources import read_text

from chimera.utils import parse_fasta
from chimera.core.domain.hmmer import HmmerDomainFinder

ctcf = read_text('chimera.data.sample', 'ctcf.fa')
ctcf = parse_fasta(ctcf)[0]


class HmmrTestCase(TestCase):
    def setUp(self):
        self.domain_finder = HmmerDomainFinder()

    def tearDown(self):
        pass

    def testHmmDomains(self):
        results = self.domain_finder.find_domains([ctcf])
        self.assertEqual(7, len(results))

        top_record = results[0]
        self.assertEqual('zf-C2H2', top_record['alihmmname'])
        # Note - Floating point values are not very accurate since we're essentially parsing Hmmer --domtblout results
        # The corresponding web api call gives us 18.923152923584
        self.assertAlmostEqual(18.9, float(top_record['bitscore']))
        self.assertEqual(266, int(top_record['alisqfrom']))
        self.assertEqual(288, int(top_record['alisqto']))
