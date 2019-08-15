from unittest import TestCase

from chimera.core.hmmr import find_hmmr_domains_web
from chimera.data.sample import ctcf


class HmmrWebTestCase(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testHmmScan(self):
        results = find_hmmr_domains_web(ctcf)
        self.assertEqual(5, len(results))

        top_domain = results[0]
        self.assertEqual('zf-H2C2_2', top_domain['name'])
        self.assertEqual(11, top_domain['ndom'])

        top_occurance = top_domain['domains'][0]
        self.assertAlmostEqual(5.89632320404053, top_occurance['bitscore'])
        self.assertEqual(262, top_occurance['alisqfrom'])
        self.assertEqual(274, top_occurance['alisqto'])
