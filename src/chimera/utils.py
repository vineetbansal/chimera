import importlib.resources
import pandas as pd

import chimera.data

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
    vdw_df = vdw_df.replace({pd.np.nan: None})


def vdw_radius(element):
    """
    Return VDW radius of element
    :param element: Element symbol (case-sensitive)
    :return: Float indicating the VDW radius
    """
    s = vdw_df.loc[element]
    return s.vdw_radius or s.ionic_radius or 1.5
