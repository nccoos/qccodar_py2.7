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

import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

from .qcutils import *
from .codarutils import *

__version__ = get_distribution("qccodar").version

def manual(datadir, pattern):
    """ Manual mode runs qc and merge on all files in datadir """

    # get file listing of datadir
    fns = recursive_glob(os.path.join(datadir, 'RadialMetric', pattern), 'RDL*.ruv')

    # handle if no files to process
    if not fns:
        print "Warn: qccodar manual --datadir %s --pattern %s" % (datadir, pattern)
        print "No files RDL*.ruv found in %s" % fulldatadir
        return
    
    print 'qccodar (manual) -- QC Processing RadialMetric to RadialShorts_qcd: ...'

    # do qc for each file in the datadir --> output to RadialShorts_qcd
    for fullfn in fns:
        print fullfn
        fn = os.path.basename(fullfn)
        do_qc(datadir, fn, pattern)

    # get file list of RadialShorts
    # depending on system and desired time span for merge, change the target time for file search
    fns = recursive_glob(os.path.join(datadir, 'RadialShorts_qcd', pattern), 'RDL*00.ruv')

    print 'qccodar (manual) -- Merging RadialShorts_qcd to Radials_qcd: ...'

    # run LLUVMerger for each
    for fullfn in fns:
        print fullfn
        fn = os.path.basename(fullfn)
        run_LLUVMerger(datadir, fn, pattern)

def auto(datadir, pattern, fullfn):
    """ Auto mode runs qc and merge when new files generated in path being watched """

    numfiles = 3
    
    # get file listing of RadialMetric folder in datadir
    indir = os.path.join(datadir, 'RadialMetric', pattern)
    fns = recursive_glob(indir, 'RDL*.ruv')

    if numfiles == 1 and len(fns)>=1:
        fullfn = fns[-1]
    elif numfiles == 3 and len(fns)>=2:
        fullfn = fns[-2]
    elif numfiles == 5 and len(fns)>=3:
        fullfn = fns[-3]
    else:
        print "qccodar (auto): Not enough data in %s to run qc" % indir
        return

    print fullfn
    fn = os.path.basename(fullfn)
    rsdfn = do_qc(datadir, fn, pattern)

    # get file listing of RadialShorts_qcd folder in datadir
    indir = os.path.join(datadir, 'RadialShorts_qcd', pattern)
    fns = recursive_glob(indir, 'RDL*00.ruv')

    if len(fns) >= 1:
        fullfn = fns[-1]
        if rsdfn and rsdfn == fullfn:
            fn = os.path.basename(fullfn)
            run_LLUVMerger(datadir, fn, pattern)

# Watcher() and Handler() classes based on
# https://www.michaelcho.me/article/using-pythons-watchdog-to-monitor-changes-to-a-directory
# and http://ginstrom.com/scribbles/2012/05/10/continuous-integration-in-python-using-watchdog/
#
# Perform any actions, using any triggers that you choose.
#

class Watcher:
    """ Watch for changes in directory and trigger an event handler """

    def __init__(self):
        self.observer = Observer()

    def run(self, datadir, pattern):
        event_handler = Handler(datadir, pattern)
        fulldatadir = os.path.join(datadir, 'RadialMetric', pattern)
        self.observer.schedule(event_handler, fulldatadir, recursive=True)
        self.observer.start()
        try:
            while True:
                time.sleep(5)
        except:
            self.observer.stop()
            print '... qccodar auto-mode stopped'

        self.observer.join()


class Handler(FileSystemEventHandler):
    """
    Specifically, handle doing qc and merge when a file is created in watched directory
    """

    def __init__(self, datadir, pattern):
        self.datadir = datadir
        self.pattern = pattern

    def on_any_event(self, event):
        if event.is_directory:
            return None
        elif event.event_type == 'created':
            # Take any action here when a file is first created.
            print "Received created event - %s." % event.src_path
            print 'datadir = %s' % self.datadir
            print 'pattern = %s' % self.pattern
            auto(self.datadir, self.pattern, event.src_path)
            
        elif event.event_type == 'modified':
            # Taken any action here when a file is modified.
            print "Received modified event - %s." % event.src_path


def main():
    """Run qccodar from the command line."""
    from docopt import docopt

    arguments = docopt(__doc__, version="qccodar %s" % __version__)
    # print arguments

    datadir, pattern = arguments['--datadir'], arguments['--pattern']
    if arguments['manual']:
        runarg = 'manual'
    elif arguments['auto']:
        runarg = 'auto'
    else:
        runarg = ''

    indatadir = os.path.join(datadir, 'RadialMetric', pattern)
    if not os.path.isdir(indatadir):
        print "Error: qccodar %s --datadir %s --pattern %s" % (runarg, datadir, pattern)
        print "Directory does not exist: %s " % indatadir
        return

    outdir1 = os.path.join(datadir, 'RadialShorts_qcd', pattern)
    if not os.path.exists(outdir1):
        os.makedirs(outdir1)

    outdir2 = os.path.join(datadir, 'Radials_qcd', pattern)
    if not os.path.isdir(outdir2):
        os.makedirs(outdir2)
 
    if arguments['manual']:
        manual(datadir, pattern)
        return

    # create watchdog to monitor datadir
    if arguments['auto']:
        w = Watcher()
        print 'Starting qccodar auto-mode (Press cntrl-C to exit)...'
        w.run(datadir, pattern)
        return

if __name__ == "__main__":
    main()
