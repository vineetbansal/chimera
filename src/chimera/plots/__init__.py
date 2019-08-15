import json
import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go

from chimera import config


def binding_freq_plot_data_domain(pfam_id):

    from chimera import df_dl

    df_dl = df_dl[
        (df_dl.num_nonidentical_instances >= config.web.min_instances) &
        (df_dl.num_structures >= config.web.min_structures) &
        (df_dl.pfam_id == pfam_id)
    ]

    bars = []
    for ligand_type in ('ion', 'metabolite', 'sm'):
        df = df_dl[df_dl.ligand_type == ligand_type]
        if len(df) == 1:
            row = df.iloc[0]
            bars.append(
                go.Bar(
                    x=list(range(1, row.domain_length + 1)),
                    y=list(map(float, row.binding_frequencies.split(','))),
                    name=ligand_type
                )
            )

    return json.dumps(bars, cls=plotly.utils.PlotlyJSONEncoder)


def binding_freq_plot_data_sequence(seq, df):

    sequence_length = len(seq)
    ligand_types = df.ligand_type.unique()

    data = np.zeros((len(ligand_types), sequence_length))
    for j, ligand_type in enumerate(ligand_types):
        bf = df[df.ligand_type == ligand_type].groupby('seq_i')['binding_frequency'].max()
        bf = bf.reindex(pd.RangeIndex(1, sequence_length + 1), fill_value=0)
        data[j, :] = bf.values

    bars = []
    for ligand_type, row in zip(ligand_types, data):
        bars.append(
            go.Bar(
                x=list(range(1, sequence_length + 1)),
                y=row,
                name=ligand_type
            )
        )

    return bars
