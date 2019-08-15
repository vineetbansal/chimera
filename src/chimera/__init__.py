__version__ = "0.1.10"


import os
import logging.config
import logging
from types import SimpleNamespace
import json
import pandas as pd
from importlib.resources import read_text, path

import chimera
import chimera.data


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
    try:
        config = json.loads(
            s,
            object_hook=lambda d: SimpleNamespace(
                **{k: (v.format(**os.environ) if isinstance(v, str) else v) for k, v in d.items()}
            )
        )
    except KeyError as e:
        raise RuntimeError(f'Environment variable {e.args[0]} missing')

    return config


config = setup_config()

# TODO: Why do we have 2 tsvs?
# TODO: Names of these dataframe directly ported from R - not intuitive!

df_bp = None
with path(chimera.data, 'interacdome_fordownload.tsv') as p:
    df_bp = pd.read_csv(p, sep='\t', header=0)
    logger.info(f'Read {len(df_bp)} records from interacdome_fordownload.tsv')

df_dl = None
with path(chimera.data, 'interacdome_allresults.tsv') as p:
    df_dl = pd.read_csv(p, sep='\t', header=0)
    logger.info(f'Read {len(df_dl)} records from interacdome_allresults.tsv')

# Unpivoted binding-frequencies table which has been pre-filtered onwedsf
# num_nonidentical_instances/num_structures/max_achieved_precision as per the config file
# TODO: Generate/cache on demand rather than ahead of time!
binding_frequencies_interacdome = None
with path(chimera.data, 'binding_frequencies_interacdome.csv') as p:
    binding_frequencies_interacdome = pd.read_csv(p, index_col=False)
    logger.info(f'Read {len(binding_frequencies_interacdome)} records from binding_frequencies_interacdome.csv')

binding_frequencies_dsprint = None
with path(chimera.data, 'binding_frequencies_dsprint.csv') as p:
    binding_frequencies_dsprint = pd.read_csv(p, index_col=False)
    logger.info(f'Read {len(binding_frequencies_dsprint)} records from binding_frequencies_dsprint.csv')
