from .version import version as __version__


import os
import logging.config
import logging
from types import SimpleNamespace
import json
import pandas as pd
from string import Formatter
from collections import defaultdict
from importlib.resources import read_text, path

import chimera
import chimera.data

from chimera.exceptions import handle_exception
import sys
sys.excepthook = handle_exception

logger = logging.getLogger(__name__)


def setup_config():
    s = read_text(chimera, 'config.json')
    d = json.loads(s)

    # logging.config doesn't support configuration from an object, but does support dictConfig,
    # so use the dict obtained from json.
    if 'logging' in d:
        logging.config.dictConfig(d['logging'])
    else:
        logging.basicConfig(level=logging.INFO)

    # Now that logging is configured, reload the json, but now with an object hook
    # so we have cleaner access to keys by way of attributes.
    # For string-type values, we try to expand them using environment variables
    # This allows us to have format-specifier style values in the json, like:
    # ..
    # "DB_URI": "mysql://{DBUSER}:{DBPASS}",
    # ..
    # For non-string type values, (this include SimpleNamespace type itself), we simply get the value as is,
    # allowing us chained attribute access).
    # This implementation uses a defaultdict to resolve these format string identifiers, and returns a blank string ''
    # for any environment variables that are missing

    default_dict = defaultdict(str)
    for k, v in os.environ.items():
        default_dict[k] = v
    fmtr = Formatter()

    config = json.loads(
        s,
        object_hook=lambda d: SimpleNamespace(
            **{k: (fmtr.vformat(v, (), default_dict) if isinstance(v, str) else v) for k, v in d.items()}
        )
    )

    return config


config = setup_config()

# The 10 different ligand-types we support for ligand-protein binding frequency algorithms
LIGAND_TYPES = ('peptide', 'ion', 'metabolite', 'sm', 'dna', 'dnabase', 'dnabackbone', 'rna', 'rnabase', 'rnabackbone')

# TODO: Names of these DataFrames directly ported from R - not intuitive!

df_dl = None
interacdome_pfam_ids = None
with path(chimera.data, 'interacdome_fordownload.tsv') as p:
    df = pd.read_csv(p, sep='\t', header=0)
    logger.info(f'Read {len(df)} records from interacdome_fordownload.tsv')

    # TODO: The way results are filtered from the 'master' tsv file for InteracDome is to compare the 'ligand_type'
    # column values with the uppercased ligand types, with a '_' appended at the end.
    # This leaves out 7 records that don't have a '_' at the end
    # (i.e. records that have ligand types 'SM', 'RNA', 'DNA')
    # We do the same here for backward compatibility, but this will need further investigation.
    old_ligand_types = [x.upper() + '_' for x in LIGAND_TYPES]
    df_dl = df[df['ligand_type'].isin(old_ligand_types)]
    df_dl['ligand_type'] = df_dl['ligand_type'].map(lambda x: x.replace('_', '').lower())
    logger.info(f'Read {len(df_dl)} filtered binding frequency records for InteracDome after filtering for ligand type')

    # Filtered records to be displayed on the website
    df_dl_filtered = df_dl[
        (df_dl.num_nonidentical_instances >= config.web.min_instances) &
        (df_dl.num_structures >= config.web.min_structures)
    ]
    logger.info(f'Read {len(df_dl_filtered)} filtered binding frequency records for InteracDome for website display')

    interacdome_pfam_ids = pd.unique(df_dl_filtered['pfam_id'])

df_dl_dsprint = None
dsprint_pfam_ids = None
with path(chimera.data, 'dsprint_fordownload.tsv') as p:
    df_dl_dsprint = pd.read_csv(p, sep='\t', header=0)
    logger.info(f'Read {len(df_dl_dsprint)} records from dsprint_fordownload.tsv')

    dsprint_pfam_ids = pd.unique(df_dl_dsprint['pfam_id'])

# Unpivoted binding-frequencies table which has been pre-filtered on
# num_nonidentical_instances/num_structures/max_achieved_precision as per the config file
# TODO: Generate/cache on demand rather than ahead of time!
binding_frequencies_interacdome = None
with path(chimera.data, 'binding_frequencies_interacdome.parquet') as p:
    binding_frequencies_interacdome = pd.read_parquet(p)
    logger.info(f'Read {len(binding_frequencies_interacdome)} records from binding_frequencies_interacdome.parquet')

binding_frequencies_dsprint = None
with path(chimera.data, 'binding_frequencies_dsprint.parquet') as p:
    binding_frequencies_dsprint = pd.read_parquet(p)
    logger.info(f'Read {len(binding_frequencies_dsprint)} records from binding_frequencies_dsprint.parquet')
