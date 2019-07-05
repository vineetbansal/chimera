import pandas as pd
from chimera.old_interacdome.calculate_distances import create_distlist_files
from chimera.distance import ANNOTATION_COLUMNS

ANNOTATION_FILE = '/media/vineetb/t5-vineetb/biolip/processed_data/annotations/current_annotations.txt'


if __name__ == '__main__':

    annot_df = pd.read_csv(ANNOTATION_FILE, sep='\t', header=None, names=ANNOTATION_COLUMNS)

    pdb_ids = list(pd.unique(annot_df['pdb_id']))
    pdb_ids = [pdb_id for pdb_id in pdb_ids if pdb_id.startswith('2lue')]

    create_distlist_files(
        ANNOTATION_FILE,
        pdb_ids,
        receptor_pdb_dir='/media/vineetb/t5-vineetb/biolip/downloaded_data/receptor/',
        ligand_pdb_dir='/media/vineetb/t5-vineetb/biolip/downloaded_data/ligand/',
        distlist_dir='/media/vineetb/t5-vineetb/biolip/processed_data/distances/',
        annot_dir='/media/vineetb/t5-vineetb/biolip/processed_data/annotations/',
        force_overwrite=True
    )
