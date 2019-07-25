from unittest import TestCase
from chimera import df_pfam


class TsvTestCase(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testNPfamIds(self):
        """
        All distinct pfam ids in the table
        """
        self.assertEqual(4128, len(df_pfam['pfam_id'].unique()))
