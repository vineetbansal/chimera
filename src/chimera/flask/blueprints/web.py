import logging
import pandas as pd
import json
import plotly
from importlib.resources import read_text
from flask import Blueprint, request, render_template, send_file, session, abort

from chimera import config
from chimera.utils import parse_fasta
from chimera.tasks import query
from chimera.plots import binding_freq_plot_data_domain, binding_freq_plot_data_sequence

bp = Blueprint('web', __name__)
logger = logging.getLogger(__name__)


ctcf = read_text('chimera.data.sample', 'ctcf.fa')


@bp.route('/dsprint', methods=['GET', 'POST'])
def dsprint():
    from chimera import df_dl_dsprint

    pfam_ids = pd.unique(df_dl_dsprint['pfam_id'])

    selected_pfam_id = None
    data = ''
    if request.method == 'POST':
        selected_pfam_id = request.form['pfam_id']
        data = binding_freq_plot_data_domain(selected_pfam_id, algorithm='dsprint')
        data = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return render_template('dsprint.html', pfam_ids=pfam_ids, selected_pfam_id=selected_pfam_id, data=data)


@bp.route('/dpuc2')
def dpuc2():
    return render_template('dpuc2.html')


@bp.route('/domstratstats')
def domstratstats():
    return render_template('domstratstats.html')


@bp.route('/interacdome', methods=['GET', 'POST'])
def interacdome():
    from chimera import df_dl

    df_dl = df_dl[
        (df_dl.num_nonidentical_instances >= config.web.min_instances) &
        (df_dl.num_structures >= config.web.min_structures)
    ]

    pfam_ids = pd.unique(df_dl['pfam_id'])

    selected_pfam_id = None
    data = ''
    if request.method == 'POST':
        selected_pfam_id = request.form['pfam_id']
        data = binding_freq_plot_data_domain(selected_pfam_id, algorithm='interacdome')
        data = json.dumps(data, cls=plotly.utils.PlotlyJSONEncoder)

    return render_template('interacdome.html', pfam_ids=pfam_ids, selected_pfam_id=selected_pfam_id, data=data)


@bp.route('/interacdome_faq')
def interacdome_faq():
    return render_template('interacdome_faq.html')


@bp.route('/interacdome_download')
def interacdome_download():
    return render_template('interacdome_download.html')


@bp.route('/', methods=['GET', 'POST'])
def index():

    error_msg = ''
    data_plotly = []
    seq_text = ''
    domain_dataframe = None
    n_hits = 0
    algorithm = ''
    domain_algorithm = ''
    email_address = ''
    full_domains = False

    job_id = ''
    n_sequences = 0
    max_sequences_interactive = config.web.max_sequences_interactive

    if request.method == 'POST':

        seq_file_text = request.files['seqFile'].read().decode('utf8')
        if seq_file_text:
            seq_text = seq_file_text
        else:
            seq_text = request.form['seqTextArea']
        domain_algorithm = request.form['algorithm0Select']
        algorithm = request.form['algorithm1Select']
        email_address = request.form['emailaddress'].strip()
        full_domains = request.form.get('fullDomainCheck', 'off') == 'on'

        sequences = parse_fasta(seq_text)
        n_sequences = len(sequences)

        # parse_fasta seems to parse '' as a single blank sequence!
        if n_sequences == 0 or str(sequences[0].seq) == '':
            error_msg = 'No sequence(s) detected.'
        elif n_sequences > max_sequences_interactive:
            if not email_address:
                error_msg = 'Please provide an email address for sending results.'
            else:
                job = query.delay(seq_text=seq_text, save_results=True, email_address=email_address,
                                  algorithm=algorithm, domain_algorithm=domain_algorithm, silent=True,
                                  full_domains=full_domains)
                job_id = job.id
        else:
            domain_dataframe, binding_dataframe, result_filename = query(sequences=sequences, save_results=True,
                                                                         email_address=email_address,
                                                                         algorithm=algorithm,
                                                                         domain_algorithm=domain_algorithm,
                                                                         full_domains=full_domains)
            session['result_filename'] = result_filename

            for i, sequence in enumerate(sequences):
                d = binding_freq_plot_data_sequence(str(sequence.seq), domain_dataframe, binding_dataframe)
                data_plotly.append(d)

            n_hits = len(domain_dataframe)

    data_plotly = json.dumps(data_plotly, cls=plotly.utils.PlotlyJSONEncoder)

    return render_template(
        'index.html',
        n_hits=n_hits,

        # For displaying warnings/messages
        error_msg=error_msg,
        job_id=job_id,
        n_sequences=n_sequences,
        max_sequences_interactive=max_sequences_interactive,

        df=domain_dataframe,
        data_plotly=data_plotly,
        seq=seq_text,
        sample_seq=ctcf,
        algorithm=algorithm,
        domain_algorithm=domain_algorithm,
        email_address=email_address,
        full_domains=full_domains
    )


@bp.route('/seq_results')
def seq_results():
    try:
        return send_file(session['result_filename'], as_attachment=True, attachment_filename='results.csv')
    except:  # nopep8
        abort(404)


