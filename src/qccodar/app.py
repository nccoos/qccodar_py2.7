"""Quality control CODAR (qccodar) RadialMetric data. 

Usage:
  qccodar (auto | manual) [options]
  qccodar (--help | --version)

Options:
  -h --help             Show this help message and exit
  --version             Show version
  -v --verbose          Verbocity
  -d DIR --datadir DIR  Data directory to process [default: /Codar/SeaSonde/Data]
  -p PAT --pattern PAT  Pattern type [default: IdealPattern]
  
"""

import os
import re
import glob

from pkg_resources import get_distribution

from .qcutils import *
from .codarutils import *

__version__ = get_distribution("qccodar").version


def manual(datadir, pattern):
    """ Manual mode runs qc on all files in datadir """

    fulldatadir = os.path.join(datadir, 'RadialMetric', pattern)

    if not os.path.isdir(fulldatadir):
        print "Error: qccodar manual --datadir %s --pattern %s" % (datadir, pattern)
        print "Directory does not exist: %s " % fulldatadir
        return
    
    # get file listing of datadir
    fns = recursive_glob(os.path.join(datadir, 'RadialMetric', pattern), 'RDL*.ruv')

    # handle if no files to process
    if not fns:
        print "Warn: qccodar manual --datadir %s --pattern %s" % (datadir, pattern)
        print "No files RDL*.ruv found in %s" % fulldatadir
        return
    
    print 'QC Processing RadialMetric to RadialShorts_qcd: ...'

    # do qc for each file in the datadir --> output to RadialShorts_qcd
    for fullfn in fns:
        print fullfn
        fn = os.path.basename(fullfn)
        do_qc(datadir, fn, pattern)

    # get file list of RadialShorts
    fns = recursive_glob(os.path.join(datadir, 'RadialShorts_qcd', pattern), 'RDL*00.ruv')

    print 'Merging RadialShorts_qcd to Radials_qcd: ...'
    outdir = os.path.join(datadir, 'Radials_qcd', pattern)

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # run LLUVMerger for each
    for fullfn in fns:
        print fullfn
        run_LLUVMerger(fullfn, outdir, pattern)


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
