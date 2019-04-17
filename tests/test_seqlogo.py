from unittest import TestCase
import importlib.resources
import pandas as pd
import seqlogo
import chimera.data.pfms


class LogoTestCase(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testSeqLogo(self):
        with importlib.resources.path(chimera.data.pfms, 'PF00013_KH_1.pfm') as path:
            df = pd.read_csv(path, sep='\t', header=None, index_col=0)
            pfm = seqlogo.Pfm(df, alphabet_type='AA')
            result = seqlogo.seqlogo(pfm, ic_scale=False, format='png', filename='logo.png')

            # As long as no exceptions were raised, we're good
            self.assertIsNone(result)
