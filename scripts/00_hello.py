import importlib.resources
import pandas as pd
import chimera.data.pfms
import seqlogo
from seqlogo import utils


if __name__ == '__main__':
    with importlib.resources.path(chimera.data.pfms, 'PF00013_KH_1.pfm') as path:
        df = pd.read_csv(path, sep='\t', header=None, index_col=0)

        # If we know that df.index.values will always be 20 and in order 'ACDEFGHIKLMNPQRSTVWY',
        # then the following 2 lines can simplified to:
        #   pfm = seqlogo.Pfm(df, alphabet_type='AA')
        #   seqlogo.seqlogo(pfm, ic_scale=False, size='xlarge', format='png', filename='logo.png', stacks_per_line=pfm.length)
        pfm = seqlogo.Pfm(df, alphabet_type='custom', alphabet=df.index.values, background=utils._AA_background)
        x = seqlogo.seqlogo(pfm, ic_scale=False, size='xlarge', format='png', filename='logo.png', stacks_per_line=pfm.length, color_scheme='hydrophobicity')
        print(x)
