from unittest import TestCase
from importlib.resources import path
from subprocess import run

import chimera.data.sample as sample
from chimera.utils import find_hmmr_domains_local
from chimera.data.sample import ctcf


class HmmrTestCase(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testHmmSearch(self):
        with path(sample, 'globins4.hmm') as hmm_path:
            with path(sample, 'globins45.fa') as fa_path:
                p = run(['hmmsearch', hmm_path, fa_path], capture_output=True)
                self.assertEquals(0, p.returncode)

    def _testHmmScan(self):
        results = find_hmmr_domains_local(ctcf)
        self.assertEqual(5, len(results))

        top_domain = results[0]
        self.assertEqual('zf-H2C2_2', top_domain['name'])
        self.assertEqual(11, top_domain['ndom'])

        top_occurance = top_domain['domains'][0]
        self.assertAlmostEqual(5.89632320404053, top_occurance['bitscore'])
        self.assertEqual(262, top_occurance['alisqfrom'])
        self.assertEqual(274, top_occurance['alisqto'])
