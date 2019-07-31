from importlib.resources import path
from unittest import TestCase
import chimera.data.sample
import chimera.data.sample.ligand
import chimera.data.sample.receptor
from chimera.distance import create_distance_file


class CalcDistanceTestCase(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testDistanceFile(self):

        with path(chimera.data.sample.ligand, '2lue_III_B_1.pdb') as ligand_filepath:
            with path(chimera.data.sample.receptor, '2lueA.pdb') as receptor_filepath:

                create_distance_file(
                    pdb_id='2lue',
                    pdb_chains=['A'],
                    receptor_filepaths=[receptor_filepath],
                    ligand_ids=['III'],
                    ligand_filepaths=[ligand_filepath],
                    distance_filepath='distance.txt',
                    include_backbone=False,
                    distance_cutoff=20,
                    compressed=False,
                    calculate_overlap=True
                )
