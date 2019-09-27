from unittest import TestCase
from importlib.resources import read_text

from chimera.utils import parse_fasta
from chimera.core.domain.domstratstats import DomStratStatsDomainFinder

ctcf = read_text('chimera.data.sample', 'ctcf.fa')
ctcf = parse_fasta(ctcf)[0]


class HmmrTestCase(TestCase):
    def setUp(self):
        self.domain_finder = DomStratStatsDomainFinder()

    def tearDown(self):
        pass

    def testHmmDomains(self):
        results = self.domain_finder.find_domains([ctcf])
        self.assertEqual(19, len(results))

