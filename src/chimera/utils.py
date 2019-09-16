import subprocess
import os.path
import logging
import importlib.resources
import pandas as pd
from io import StringIO
from Bio import SeqIO

import smtplib
import ssl
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

import chimera.data
from chimera import config

logger = logging.getLogger(__name__)

vdw_df = None
with importlib.resources.path(chimera.data, 'vdw.txt') as path:
    vdw_df = pd.read_csv(
        path,
        sep='\t',
        header=None,
        comment='#',
        names=(
            'atomic_number',
            'element_symbol',
            'atomic_radius',
            'ionic_radius',
            'covalent_radius',
            'vdw_radius',
            'crystal_radius'
        ),
        index_col=1
    )
    vdw_df.index = vdw_df.index.str.lower()
    vdw_df = vdw_df.replace({pd.np.nan: None})


def vdw_radius(element):
    """
    Return VDW radius of element
    :param element: Element symbol (case-sensitive)
    :return: Float indicating the VDW radius
    """
    s = vdw_df.loc[element.lower()]
    return s.vdw_radius or s.ionic_radius or 1.5


_ligand_to_groups = {}


def ligand_to_groups(ligand_id):

    global _ligand_to_groups
    if not _ligand_to_groups:
        with importlib.resources.path(chimera.data, 'ligand_groups.txt') as path:
            df = pd.read_csv(
                path,
                sep='\t',
                header=None,
                comment='#',
                names=('group_name', 'ligand_id'),
                na_filter=False  # important since 'NA' is a valid ligand identifier!
            )

            _ligand_to_groups = {}
            for l_id, _df in df.groupby('ligand_id'):
                group_names = set(_df.group_name.values)

                # a single molecule must belong exclusively to the nucleic acid, ion, or small molecule groups
                if 'NUCACID_' in group_names:
                    [group_names.discard(g) for g in ['ION_', 'METABOLITE_', 'DRUGLIKE_', 'SM_']]
                elif 'ION_' in group_names:
                    [group_names.discard(g) for g in ['METABOLITE_', 'DRUGLIKE_', 'SM_']]

                _ligand_to_groups[l_id] = sorted(group_names)

    group_names = ['ALL_', ligand_id] + _ligand_to_groups[ligand_id]
    if not ('NUCACID_' in group_names or 'ION_' in group_names or 'III' in group_names):
        group_names.append('SM_')

    # Last-minute Ligand Group renamings
    group_renamings = {
        'NUCDNA': 'DNABASE_',
        'NUCDNAB': 'DNABACKBONE_',
        'NUCRNA': 'RNABASE_',
        'NUCRNAB': 'RNABACKBONE_',
        'III': 'PEPTIDE_'
    }

    group_names = [group_renamings.get(g, g) for g in group_names]

    return group_names


def parse_fasta(s, default_seq_id='seq0'):
    sequences = list(SeqIO.parse(StringIO(s), 'fasta'))
    if not sequences:
        sequences = list(SeqIO.parse(StringIO(f'>{default_seq_id}\n' + s), 'fasta'))
    return sequences


def get_full_version():
    """
    Get as much version information as we can, including git info (if applicable)
    This method should never raise exceptions!
    :return: A version no. in the form:
        <maj>.<min>.<bld>
            If we're running as a package distributed through setuptools
        <maj>.<min>.<bld>.<rev>
            If we're running as a 'regular' python source folder, possibly locally modified
            <rev> is one of:
                'src': The package is running as a source folder
                <git_tag> or <git_rev> or <git_rev>-dirty: A git tag or commit revision, possibly followed by a suffix
                    '-dirty' if source is modified locally
                'x':   The revision cannot be determined
    """
    import chimera
    from chimera.version import version

    rev = ''
    try:
        path = chimera.__path__[0]
        if os.path.isdir(path):
            # We have a package folder where we can get git information
            try:
                rev = subprocess.check_output(['git', 'describe', '--tags', '--always', '--dirty'], stderr=subprocess.STDOUT, cwd=path).decode('utf-8').strip()
            except (FileNotFoundError, subprocess.CalledProcessError):
                # no git or not a git repo? assume 'src'
                rev = 'src'
        else:
            # We're very likely running as a package
            rev = ''
    except:
        # Something unexpected happened - rev no. defaults to 'x'
        rev = 'x'

    if rev == '':
        return version
    else:
        return f'{version}.{rev}'


def email(to, subject, body, attachments=None):

    message = MIMEMultipart()
    message["From"] = config.mail.sender
    message["To"] = to
    message["Subject"] = subject

    message.attach(MIMEText(body, "plain"))

    if attachments is not None:
        for attachment in attachments:
            with open(attachment, "rb") as f:
                part = MIMEBase("application", "octet-stream")
                part.set_payload(f.read())
            encoders.encode_base64(part)

            part.add_header(
                "Content-Disposition",
                f"attachment; filename= {os.path.basename(attachment)}",
            )
            message.attach(part)

    text = message.as_string()

    context = ssl.create_default_context()
    with smtplib.SMTP_SSL(config.mail.smtp_server, int(config.mail.smtp_port), context=context) as server:
        server.login(config.mail.smtp_username, config.mail.smtp_password)
        server.sendmail(config.mail.sender, to, text)
