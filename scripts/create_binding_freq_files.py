import pandas as pd
from chimera import df_dl, config


if __name__ == '__main__':

    df = df_dl

    df = df[
        (df['num_nonidentical_instances'] >= config.web.min_instances) &
        (df['num_structures'] >= config.web.min_structures) &
        (df['max_achieved_precision'] >= config.web.min_achieved_precision)
    ]

    l = []
    for pfam_id, _df in df.groupby('pfam_id'):
        for _, row in _df.iterrows():
            for i, bf in enumerate(map(float, row.binding_frequencies.split(',')), start=1):
                l.append({
                    'pfam_id': pfam_id,
                    'match_state': i,
                    'ligand_type': row.ligand_type,
                    'binding_frequency': bf
                })

    pfam_df = pd.DataFrame(l)
    pfam_df.to_csv('binding_frequencies.csv', index=False)
