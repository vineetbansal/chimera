from unittest import TestCase
import pandas as pd
from chimera import df_dl


class TsvTestCase(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testNDomains(self):
        """
        No. of unique domain identifier rows that we have
        """
        self.assertEqual(4128, len(df_dl['pfam_id'].unique()))

    def testDomainLength(self):
        """
        Each unique domain identifier has a unique domain length
        """
        for pfam_id, _df in df_dl.groupby('pfam_id'):
            self.assertEqual(1, len(_df['domain_length'].unique()))

    def testNVisibleDomains(self):
        """
        All unique domain identifiers that are displayed on the website dropdown
        """
        df = df_dl[(df_dl.num_nonidentical_instances >= 3) & (df_dl.num_structures >= 3)]
        pfam_ids = pd.unique(df['pfam_id'])
        self.assertEqual(2263, len(pfam_ids))

    def testNBindingFreqs(self):
        """
        No. of binding frequencies in each row is equal to the domain length
        :return:
        """
        for i, row in df_dl.iterrows():
            self.assertEqual(row.domain_length, len(row.binding_frequencies.split(',')))

