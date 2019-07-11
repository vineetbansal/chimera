import os
import gzip
import pandas as pd

from chimera.distance import ANNOTATION_COLUMNS, create_distance_file

ANNOTATION_FILE = '/media/vineetb/t5-vineetb/biolip/processed_data/annotations/current_annotations.txt'
COMPRESS = False

if __name__ == '__main__':
    receptor_dir = '/media/vineetb/t5-vineetb/biolip/downloaded_data/receptor/'
    ligand_dir = '/media/vineetb/t5-vineetb/biolip/downloaded_data/ligand/'

    annot_df = pd.read_csv(ANNOTATION_FILE, sep='\t', header=None, names=ANNOTATION_COLUMNS)

    pdb_ids = annot_df.pdb_id.unique().tolist()
    # TEMPORARY
    pdb_ids = [pdb_id for pdb_id in pdb_ids if pdb_id.startswith('2l')]

    for pdb_id in pdb_ids:
        print(f'Processing PDB {pdb_id}')
        df = annot_df[annot_df.pdb_id == pdb_id]

        pdb_chains = df.pdb_chain.tolist()
        ligand_ids = df.ligand_id.tolist()
        ligand_chains = df.ligand_chain.tolist()
        ligand_snos = df.ligand_serial_number.tolist()

        try:
            create_distance_file(
                pdb_id=pdb_id,
                pdb_chains=pdb_chains,
                receptor_filepaths=[os.path.join(receptor_dir, f'{pdb_id}{pdb_chain}.pdb') for pdb_chain in pdb_chains],
                ligand_ids=ligand_ids,
                ligand_filepaths=[os.path.join(ligand_dir, f'{pdb_id}_{ligand_id}_{ligand_chain}_{ligand_sno}.pdb') for ligand_id, ligand_chain, ligand_sno in zip(ligand_ids, ligand_chains, ligand_snos)],
                distance_filepath='distance.txt',
                include_backbone=False,
                distance_cutoff=20,
                calculate_overlap=True,
                compressed=COMPRESS
            )

            if COMPRESS:
                s = gzip.open('distance.txt.gz', 'rb').read().decode('utf8')
            else:
                s = open('distance.txt', 'rb').read().decode('utf8')

            s1 = gzip.open(os.path.join('/media/vineetb/t5-vineetb/biolip/processed_data/distances_overlap/', pdb_id[0], pdb_id[:2], f'{pdb_id}_distances.txt.gz'), 'rb').read().decode('utf8')
            if (s==s1):
                print(pdb_id + ' GOOD')
            else:
                print(pdb_id + ' BAD')
        except AssertionError as e:
            print(pdb_id + ' BAD ' + str(e))
