from unittest import TestCase
from importlib.resources import path
from subprocess import run
import chimera.data.sample as sample


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
