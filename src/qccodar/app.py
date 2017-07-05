"""Quality control CODAR (qccodar) RadialMetric data. 

Usage:
  qccodar (auto | manual) [options]
  qccodar (--help | --version)

Options:
  -h --help             Show this help message and exit
  --version             Show version
  -d DIR --datadir DIR  Data directory to process [default: /Codar/SeaSonde/Data]
  -p PAT --pattern PAT  Pattern type [default: IdealPattern]

  
"""

import os
import re

from pkg_resources import get_distribution

from .qcutils import *
from .codarutils import *

__version__ = get_distribution("qccodar").version


def manual(datadir, pattern):
    """ Manual mode runs qc on all files in datadir """
    
    # get file listing of datadir
    fns = recursive_glob(os.path.join(datadir, 'RadialMetric', pattern), 'RDL*.ruv')
    print 'Processing: ...'

    # do qc for each file in the datadir
    for fullfn in fns:
        print fullfn
        fn = os.path.basename(fullfn)
        do_qc(datadir, fn, pattern)


def main():
    """Run qccodar from the command line."""
    from docopt import docopt

    arguments = docopt(__doc__, version="qccodar %s" % __version__)
    # print arguments

    datadir, pattern = arguments['--datadir'], arguments['--pattern']

    if arguments['manual']:
        manual(datadir, pattern)
        return

    # create watchdog to monitor datadir

if __name__ == "__main__":
    main()
