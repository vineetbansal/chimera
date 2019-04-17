from unittest import TestCase
import importlib.resources
import pandas as pd
import chimera.data


class TsvTestCase(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testNDomains(self):
        """
        All domain identifier rows that we have
        """
        with importlib.resources.path(chimera.data, 'interacdome_allresults.tsv') as path:
            df = pd.read_csv(path, sep='\t', header=0)
            self.assertEqual(11655, len(df))

    def testNVisibleDomains(self):
        """
        All domain identifiers that are displayed on the website
        """
        with importlib.resources.path(chimera.data, 'interacdome_allresults.tsv') as path:
            df = pd.read_csv(path, sep='\t', header=0)
            df = df[(df.num_nonidentical_instances >= 3) & (df.num_structures >= 3)]
            pfam_ids = pd.unique(df['pfam_id'])
            self.assertEqual(2263, len(pfam_ids))

