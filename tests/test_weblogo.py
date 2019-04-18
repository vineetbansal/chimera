from unittest import TestCase
import importlib.resources
import pandas as pd

from weblogo import LogoOptions, LogoData, LogoFormat, ColorScheme
from weblogo.seq import unambiguous_protein_alphabet
from weblogo.logo_formatter import png_formatter

import chimera.data.pfms


class LogoTestCase(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testSeqLogo(self):
        with importlib.resources.path(chimera.data.pfms, 'PF00013_KH_1.pfm') as path:
            df = pd.read_csv(path, sep='\t', header=None, index_col=0)

            # LogoData.from_counts expects counts as an array of seq_len x |alphabet|
            data = LogoData.from_counts(alphabet=unambiguous_protein_alphabet, counts=df.values.T)

            image = png_formatter(
                data,
                LogoFormat(
                    data,
                    LogoOptions(
                        fineprint=False,
                        logo_title='Title Here',
                        color_scheme=ColorScheme(),
                        stacks_per_line=data.length  # all data on one line
                    )
                )
            )

            # image has type 'bytes', write to file using:
            # open('img.png', 'wb').write(image)

            # As long as no exceptions were raised, we're good
            self.assertIsNotNone(image)

