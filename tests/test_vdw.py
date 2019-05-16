from unittest import TestCase
from chimera.utils import vdw_radius


class VDWTestCase(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testVDWH(self):
        # first entry
        self.assertEqual(1.2, vdw_radius('H'))

    def testVDWCl(self):
        # random entry
        self.assertEqual(1.75, vdw_radius('Cl'))

    def testVDWCa(self):
        # missing vdw_radius
        self.assertEqual(1.8, vdw_radius('Ca'))

    def testVDWFr(self):
        # missing ionic_radius and vdw_radius
        self.assertEqual(1.5, vdw_radius('Fr'))
