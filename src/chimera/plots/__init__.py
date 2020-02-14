import json
import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
from plotly.colors import DEFAULT_PLOTLY_COLORS

from chimera import config, LIGAND_TYPES


# Each ligand 'type' is denoted in it's own color (in the bars/legend), to conform to the colors used in publication.
LIGAND_TYPE_COLORS = {
    'rna': '#e41a1c',
    'dna': '#377eb8',
    'ion': '#4daf4a',
    'peptide': '#984ea3',
    'sm': '#ff7f00',
    'metabolite': '#9993d5',
    'dnabase': '#d8c12e',
    'dnabackbone': '#e0948e',
    'rnabase': '#7a9baf',
    'rnabackbone': '#8baabc'
}


def binding_freq_plot_data_domain(pfam_id, algorithm):

    if algorithm == 'interacdome':
        from chimera import df_dl
        df_dl = df_dl[
            (df_dl.num_nonidentical_instances >= config.web.min_instances) &
            (df_dl.num_structures >= config.web.min_structures)
        ]
    elif algorithm == 'dsprint':
        from chimera import df_dl_dsprint as df_dl

    df_dl = df_dl[df_dl.pfam_id == pfam_id]

    traces = []
    colorway = []
    for ligand_type in LIGAND_TYPES:
        df = df_dl[df_dl.ligand_type == ligand_type]
        if len(df) == 1:
            row = df.iloc[0]
            traces.append(
                go.Bar(
                    x=list(range(1, row.domain_length + 1)),
                    y=list(map(float, row.binding_frequencies.split(','))),
                    name=ligand_type
                )
            )
            colorway.append(LIGAND_TYPE_COLORS.get(ligand_type, '#ccc'))

    return {
        'data': traces,
        'layout': {
            "colorway": colorway
        }
    }


def binding_freq_plot_data_sequence(seq, domain_df, df):

    sequence_length = len(seq)
    ligand_types = df.ligand_type.unique()

    data = np.zeros((len(ligand_types), sequence_length))
    for j, ligand_type in enumerate(ligand_types):
        bf = df[df.ligand_type == ligand_type].groupby('seq_i')['binding_frequency'].max()
        bf = bf.reindex(pd.RangeIndex(1, sequence_length + 1), fill_value=0)
        data[j, :] = bf.values

    traces = []
    colorway = []
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
        colorway.append(LIGAND_TYPE_COLORS.get(ligand_type, '#ccc'))

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

    return {
        'data': traces,
        'layout': {
            "xaxis": {"anchor": "y"},
            "yaxis": {"anchor": "x", "domain": [0.25, 1.0]},
            "yaxis": {"anchor": "x", "domain": [0.25, 1.0]},
            "xaxis2": {"anchor": "y2", "matches": "x"},
            "yaxis2": {"anchor": "x2", "domain": [0.0, 0.2]},
            "colorway": colorway
        }
    }

