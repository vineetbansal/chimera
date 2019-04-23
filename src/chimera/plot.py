import io
import base64
import importlib.resources
import importlib.resources
import pandas as pd
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import chimera.data

with importlib.resources.path(chimera.data, 'interacdome_allresults.tsv') as path:
    RESULTS = pd.read_csv(path, sep='\t', header=0)
    RESULTS = RESULTS[
        (RESULTS.num_nonidentical_instances >= 3) & (RESULTS.num_structures >= 3)
    ]


def figure(x, y, title=None, title_fontsize=None, figsize=(50, 3), **kwargs):
    fig = Figure(figsize=figsize)

    # Artists to consider when determining the plot size
    artists = []

    if title is not None:
        title = fig.suptitle(title, fontsize=title_fontsize)
        artists.append(title)

    _ = FigureCanvas(fig)
    axis = fig.add_subplot(111)

    xlabel = kwargs.get('xlabel')
    if xlabel is not None:
        axis.set_xlabel(xlabel)

    ylabel = kwargs.get('ylabel')
    if ylabel is not None:
        axis.set_ylabel(ylabel)

    xticks = kwargs.get('xticks')
    if xticks is not None:
        axis.set_xticks(xticks)
        axis.xaxis.set_tick_params(rotation=45, labelsize=6)

    yticks = kwargs.get('yticks')
    if xticks is not None:
        axis.set_yticks(yticks)

    xlim = kwargs.get('xlim')
    if xlim is not None:
        axis.set_xlim(xlim)

    ylim = kwargs.get('ylim')
    if ylim is not None:
        axis.set_ylim(ylim)

    axis.grid(alpha=0.3)

    axis.bar(
        x,
        y
    )

    # lgd = axis.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08))
    # artists.append(lgd)

    img = io.BytesIO()
    fig.savefig(img, format='png', bbox_extra_artists=artists, bbox_inches='tight')

    return base64.b64encode(img.getvalue()).decode()


def binding_freq_figure(pfam_id, ligand_type, title=None, raise_errors=True, detailed=False):
    df = RESULTS[(RESULTS.pfam_id == pfam_id) & (RESULTS.ligand_type == ligand_type)]
    if len(df) != 1:
        if raise_errors:
            raise RuntimeError('Did not locate a single row of results')
    else:
        row = df.iloc[0]
        if title is None:
            title = '{} per-position binding frequencies'.format(ligand_type)

        if detailed:
            kwargs = dict(
                x=list(range(1, row.domain_length + 1)),
                y=list(map(float, row.binding_frequencies.split(','))),
                title=title,
                title_fontsize=10,
                xlabel='Position in Domain',
                ylabel='Binding Frequency',
                xlim=(0, row.domain_length + 1),
                ylim=(0, 1.01),
                xticks=list(np.arange(1, row.domain_length + 1, 1)),
                yticks=(0, 0.25, 0.5, 0.75, 1),
            )
        else:
            kwargs = dict(
                x=list(range(1, row.domain_length + 1)),
                y=list(map(float, row.binding_frequencies.split(','))),
                title=title,
                title_fontsize=20,
                xlabel=None,
                ylabel=None,
                xlim=(0, row.domain_length + 1),
                ylim=(0, 1.01),
                yticks=(0, 1),
            )

        return figure(**kwargs)

