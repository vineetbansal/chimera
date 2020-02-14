from unittest import TestCase
from importlib.resources import read_text

from chimera.utils import parse_fasta
from chimera.core.domain.hmmerweb import HmmerWebDomainFinder


ctcf = read_text('chimera.data.sample', 'ctcf.fa')
ctcf = parse_fasta(ctcf)[0]


class HmmrWebTestCase(TestCase):
    def setUp(self):
        self.domain_finder = HmmerWebDomainFinder()

    def tearDown(self):
        pass

    def testHmmScan(self):
        results = self.domain_finder.find_domains([ctcf])
        self.assertEqual(7, len(results))

        top_record = results[0]
        self.assertEqual('zf-C2H2', top_record['alihmmname'])
        self.assertAlmostEqual(18.923152923584, float(top_record['bitscore']))
        self.assertEqual(266, int(top_record['alisqfrom']))
        self.assertEqual(288, int(top_record['alisqto']))
