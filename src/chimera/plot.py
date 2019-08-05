import io
import base64
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost

from chimera import config


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

    from chimera import df_dl

    df_dl = df_dl[
        (df_dl.num_nonidentical_instances >= config.web.min_instances) &
        (df_dl.num_structures >= config.web.min_structures) &
        (df_dl.pfam_id == pfam_id) &
        (df_dl.ligand_type == ligand_type)
    ]

    if len(df_dl) != 1:
        if raise_errors:
            raise RuntimeError('Did not locate a single row of results')
    else:
        row = df_dl.iloc[0]
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


def sequence_binding_freq_figure(seq, ligand_labels, data, algorithm):

    M = len(ligand_labels)
    L = len(seq)

    fig = Figure(figsize=(200, M))
    _ = FigureCanvas(fig)
    ax = SubplotHost(fig, 1, 1, 1)
    fig.add_subplot(ax)

    fig.suptitle(algorithm)
    ax.imshow(data, cmap='gray_r', vmin=0, vmax=1, origin='lower')

    ax.yaxis.set_ticks(list(range(len(ligand_labels))))
    ax.set_yticklabels(ligand_labels)

    x_ticks = np.insert(np.arange(-1, L, 5)[1:], 0, [0])
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([i + 1 for i in x_ticks])

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(np.arange(0, L))
    ax2.set_xticklabels(seq)

    img = io.BytesIO()
    fig.savefig(img, format='png', bbox_extra_artists=[], bbox_inches='tight')

    return base64.b64encode(img.getvalue()).decode()
