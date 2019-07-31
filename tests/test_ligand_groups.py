from unittest import TestCase
from chimera.utils import ligand_to_groups


class LigandGroupsTestCase(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testLCO(self):
        groups = ligand_to_groups('LCO')
        # Ligand 'LCO' is both 'DRUGLIKE_' and 'ION_', but we filter out 'DRUGLIKE_' as per the directive:
        # "a single molecule must belong exclusively to the nucleic acid, ion, or small molecule groups"
        self.assertEqual(groups, ['ALL_', 'LCO', 'ION_'])

    def test14Y(self):
        groups = ligand_to_groups('14Y')
        self.assertEqual(groups, ['ALL_', '14Y', 'METABOLITE_', 'SM_'])
