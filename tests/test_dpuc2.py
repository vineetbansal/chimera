import os.path
from unittest import TestCase
from importlib.resources import read_text

from chimera.dpuc2 import Dpuc2

ctcf = read_text('chimera.data.sample', 'ctcf.fa')


class Dpuc2TestCase(TestCase):
    def setUp(self):
        self.dpuc2 = Dpuc2()

    def tearDown(self):
        pass

    def testScan(self):
        result_string = self.dpuc2.run_scan(ctcf)
        self.assertEqual(
            result_string,
            open(os.path.join(os.path.dirname(__file__), 'dpuc2_ctcf_results.txt'), 'r').read()
        )

    def testCsv(self):
        lines = open(os.path.join(os.path.dirname(__file__), 'dpuc2_ctcf_results.txt'), 'r').readlines()

        results = []  # A list-of-dicts that we'll convert to a DataFrame
        for line in lines[2:]:  # skip first 2 lines
            parts = line.split()
            print(len(parts))

            alihmmacc = parts[0]
            alihmmname = parts[1]
            pfam_domain = alihmmacc[:7] + '_' + alihmmname
            target_start = parts[17]
            target_end = parts[18]
            hmm_start = parts[15]
            hmm_end = parts[16]
            domain_length = parts[2]
            bit_score = 0
            reported = False
            e_value = parts[6]


