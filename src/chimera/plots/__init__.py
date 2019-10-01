import json
import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
from plotly.colors import DEFAULT_PLOTLY_COLORS

from chimera import config


def binding_freq_plot_data_domain(pfam_id):

    from chimera import df_dl

    df_dl = df_dl[
        (df_dl.num_nonidentical_instances >= config.web.min_instances) &
        (df_dl.num_structures >= config.web.min_structures) &
        (df_dl.pfam_id == pfam_id)
    ]

    bars = []
    for ligand_type in ('peptide', 'ion', 'metabolite', 'sm', 'dna', 'dnabase', 'dnabackbone', 'rna', 'rnabase', 'rnabackbone'):
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


def binding_freq_plot_data_sequence(seq, domain_df, df):

    sequence_length = len(seq)
    ligand_types = df.ligand_type.unique()

    data = np.zeros((len(ligand_types), sequence_length))
    for j, ligand_type in enumerate(ligand_types):
        bf = df[df.ligand_type == ligand_type].groupby('seq_i')['binding_frequency'].max()
        bf = bf.reindex(pd.RangeIndex(1, sequence_length + 1), fill_value=0)
        data[j, :] = bf.values

    traces = []
    for ligand_type, row in zip(ligand_types, data):
        traces.append(
            go.Bar(
                x=list(range(1, sequence_length + 1)),
                y=row,
                name=ligand_type,
                xaxis="x1",
                yaxis="y1"
            )
        )

    """
    Box plots to show domains on an x-axis that is aligned with the above x-axis for the bar plot.
    Note the following choices that affect the box plot behavior.
    'name': Plotly vertically aligns box plot(s) with distinct 'name' in its own row.
            We use the same 'name' so that all box plots are in the same line.
    'line_width': By making this 0, we prevent the 'median' line to be displayed inside each individual box.
    """
    for i, (pfam_domain, df) in enumerate(domain_df.groupby('pfam_domain')):
        for _, row in df.iterrows():
            traces.append(
                go.Box(
                    x=[row.target_start, row.target_end],
                    name=pfam_domain,
                    xaxis="x2",
                    yaxis="y2",
                    showlegend=False,
                    line_width=0,
                    fillcolor=DEFAULT_PLOTLY_COLORS[i % len(DEFAULT_PLOTLY_COLORS)]
                )
            )

    return traces
