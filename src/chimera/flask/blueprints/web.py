import importlib.resources
import base64
import pandas as pd
from flask import Blueprint, request, render_template

from weblogo import LogoOptions, LogoData, LogoFormat, ColorScheme
from weblogo.seq import unambiguous_protein_alphabet
from weblogo.logo_formatter import png_formatter

import chimera.data.pfms
from chimera.plot import RESULTS, binding_freq_figure

bp = Blueprint('web', __name__)


@bp.route('/', methods=['GET', 'POST'])
def index():
    pfam_ids = pd.unique(RESULTS['pfam_id'])

    imgs = []
    selected_pfam_id = None
    if request.method == 'POST':
        selected_pfam_id = request.form['pfam_id']
        imgs = images(selected_pfam_id, detailed=False)

    return render_template('index.html', pfam_ids=pfam_ids, selected_pfam_id=selected_pfam_id, imgs=imgs)


@bp.route('/detail/<pfam_id>')
def detail(pfam_id):
    return render_template('detail.html', pfam_id=pfam_id, imgs=images(pfam_id, detailed=True))


def images(pfam_id, detailed=True):
    imgs = list(filter(
        None,
        [binding_freq_figure(pfam_id, _type, raise_errors=False, detailed=detailed) for _type in ('ion', 'metabolite', 'sm')]
    ))

    with importlib.resources.path(chimera.data.pfms, pfam_id + '.pfm') as path:
        df = pd.read_csv(path, sep='\t', header=None, index_col=0)

        data = LogoData.from_counts(alphabet=unambiguous_protein_alphabet, counts=df.values.T)

        img = png_formatter(
            data,
            LogoFormat(
                data,
                LogoOptions(
                    fineprint=False,
                    logo_title='',
                    color_scheme=ColorScheme(),
                    stacks_per_line=data.length  # all data on one line
                )
            )
        )

        img = base64.b64encode(img).decode()
        imgs.append(img)
    return imgs


@bp.route('/faqs')
def faqs():
    return render_template('faqs.html')
