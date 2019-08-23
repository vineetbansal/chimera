import os.path
import logging
import pandas as pd
import json
import plotly
from tempfile import NamedTemporaryFile
from importlib.resources import read_text
from flask import Blueprint, request, render_template, send_file, session, abort

from chimera import config
from chimera.utils import parse_fasta
from chimera.core import query
from chimera.plots import binding_freq_plot_data_domain, binding_freq_plot_data_sequence

bp = Blueprint('web', __name__)
logger = logging.getLogger(__name__)


ctcf = read_text('chimera.data.sample', 'ctcf.fa')


@bp.route('/dsprint')
def dsprint():
    return render_template('dsprint.html')


@bp.route('/dpuc2')
def dpuc2():
    return render_template('dpuc2.html')


@bp.route('/interacdome', methods=['GET', 'POST'])
def interacdome():
    from chimera import df_dl

    df_dl = df_dl[
        (df_dl.num_nonidentical_instances >= config.web.min_instances) &
        (df_dl.num_structures >= config.web.min_structures)
    ]

    pfam_ids = pd.unique(df_dl['pfam_id'])

    selected_pfam_id = None
    data = []
    if request.method == 'POST':
        selected_pfam_id = request.form['pfam_id']
        data = binding_freq_plot_data_domain(selected_pfam_id)

    return render_template('interacdome.html', pfam_ids=pfam_ids, selected_pfam_id=selected_pfam_id, data=data)


@bp.route('/faqs')
def faqs():
    return render_template('faqs.html')


@bp.route('/', methods=['GET', 'POST'])
def index():

    data_plotly = []
    algorithm = ''
    seq_text = ''
    df = None
    n_hits = 0

    n_sequences = 0
    n_sequences_discarded = 0
    n_graphs = 0
    n_graphs_discarded = 0
    max_sequences = config.web.max_sequences
    max_graphs = config.web.max_graphs

    if request.method == 'POST':

        seq_file_text = request.files['seqFile'].read().decode('utf8')
        if seq_file_text:
            seq_text = seq_file_text
        else:
            seq_text = request.form['seqTextArea']
        domain_algorithm = request.form['algorithm0Select']
        algorithm = request.form['algorithm1Select']

        if domain_algorithm == 'dPUC2':
            from chimera.dpuc2 import Dpuc2
            s = Dpuc2().run_scan(seq_text)
            return "<div style=\"white-space: pre-wrap; font-family:monospace;\">" + s + "</div>"

        sequences = parse_fasta(seq_text)
        n_sequences_discarded = max(0, len(sequences)-max_sequences)
        if n_graphs_discarded > 0:
            sequences = sequences[:-n_sequences_discarded]
        n_sequences = len(sequences)
        n_graphs_discarded = max(0, n_sequences-max_graphs)
        n_graphs = n_sequences - n_graphs_discarded

        domain_dataframes = []
        binding_dataframes = []
        for i, sequence in enumerate(sequences):
            seq_name = sequence.name
            seq = str(sequence.seq)

            domain_dataframe, binding_dataframe = query(seq, algorithm)

            domain_dataframe['seq_index'] = i + 1
            domain_dataframe['seq_name'] = seq_name
            domain_dataframes.append(domain_dataframe)

            binding_dataframe['seq_index'] = i + 1
            binding_dataframe['seq_name'] = seq_name
            binding_dataframes.append(binding_dataframe)

            if i < n_graphs:
                bars = binding_freq_plot_data_sequence(seq, binding_dataframe)
                data_plotly.append({'bars': bars, 'seq_index': i + 1, 'seq_name': seq_name})

        df = pd.concat(domain_dataframes, axis=0)
        n_hits = len(df)

        # Create a temporary file for results
        with NamedTemporaryFile(delete=False) as f:
            pd.concat(binding_dataframes, axis=0).to_csv(f.name, index=False)
            session['result_filename'] = f.name

    data_plotly = json.dumps(data_plotly, cls=plotly.utils.PlotlyJSONEncoder)

    return render_template(
        'index.html',
        n_hits=n_hits,

        # For displaying warnings
        n_sequences=n_sequences,
        n_sequences_discarded=n_sequences_discarded,
        n_graphs=n_graphs,
        n_graphs_discarded=n_graphs_discarded,
        max_sequences=max_sequences,
        max_graphs=max_graphs,

        df=df,
        data_plotly=data_plotly,
        seq=seq_text,
        sample_seq=ctcf,
        algorithm=algorithm
    )


@bp.route('/seq_results')
def seq_results():
    try:
        return send_file(session['result_filename'], as_attachment=True, attachment_filename='results.csv')
    except:  # nopep8
        abort(404)


@bp.route('/data')
def data():
    s = ''
    data_path = request.args.get('path', '/data')
    s += '<br/>Data Path=' + data_path
    try:
        files = os.listdir(data_path)
        s += '<br/>No of files=' + str(len(files))
        for file in files:
            s += '<br/>' + str(file)
    except Exception as e:
        s += '<br/>' + str(e)

    try:
        x = open(os.path.join(data_path, 'greeting.txt'), 'r').read()
        s += '<br/>Greeting=' + str(x)
    except Exception as e:
        s += '<br/>' + str(e)

    return s
