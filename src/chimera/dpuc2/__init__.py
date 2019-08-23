import os.path
import logging
from subprocess import run
from tempfile import NamedTemporaryFile

from chimera import config

logger = logging.getLogger(__name__)


class Dpuc2:

    def __init__(self):
        self.pfam_hmm_path = os.path.join(config.dirs.data.pfam, 'Pfam-A.hmm')
        self.pfam_hmm_dat_path = os.path.join(config.dirs.data.pfam, 'Pfam-A.hmm.dat')
        self.hmmscan_bin = os.path.join(config.dirs.bin.hmmr, 'hmmscan')

    def run_scan(self, sequence):

        f_in = NamedTemporaryFile(delete=False)
        f_intermediate = NamedTemporaryFile(delete=False)
        f_out = NamedTemporaryFile(delete=False)

        f_in.write(sequence.encode('utf8'))
        f_in.close()
        logger.info(f'Wrote sequence to temporary file = {f_in.name}')

        try:
            p1 = run([
                'perl',
                '0runHmmscan.pl',
                self.hmmscan_bin,
                self.pfam_hmm_path,
                f_in.name,
                f_intermediate.name
            ], cwd=config.dirs.bin.dpuc2, capture_output=True)

            if p1.returncode != 0:
                raise RuntimeError('Perl execution failed')

            logger.info(f'Wrote dpuc2 hmmscan results to temporary file = {f_intermediate.name}')

            p2 = run([
                'perl',
                '1dpuc2.pl',
                self.pfam_hmm_dat_path,
                config.files.data.dpuc2_net,
                f_intermediate.name,
                f_out.name,
                '--noGzip'
            ], cwd=config.dirs.bin.dpuc2, capture_output=True)

            if p2.returncode != 0:
                raise RuntimeError('Perl execution failed')

            logger.info(f'Wrote dpuc2 final results to temporary file = {f_out.name}')

            return open(f_out.name).read()

        finally:
            os.unlink(f_in.name)
            os.unlink(f_intermediate.name)
            os.unlink(f_out.name)



