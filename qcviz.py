#!/usr/bin/env python
#
# Last modified: Time-stamp: <2015-06-05 17:36:23 Sara>
""" Vizualization tool for QC process and settings

For a given bearing:
  First plot -- velocities with range cell for a given bearing
  Second plot -- SNR with range or other selected radialmetric parameter
  Third plot -- compass plot of current bearing direction 
Animate with bearing.

Other GUI:
  Change threshold values for each threhold test [5.0 50.0 5.0]
  Change numfiles 3 (odd int)
  Change numdegrees 3 (odd int)
  Select different weighting factor MP (MP, SNR, NONE)

"""
from qcutils import *
from codarutils import *

import matplotlib 
print 'matplotlib backend: %s' % (matplotlib.get_backend(),)
import matplotlib.pyplot as plt

# pylab.rcParams['figure.figsize']= (11.0, 8.0)
fig, axs = plt.subplots(3,1)

axs[0].axhline(y=0, linewidth=1, color='k')
axs[0].set_ylim(-150, 150)
axs[0].set_xlabel('Range Cell')
axs[0].set_ylabel('Radial Velocity (cm/s)')
ld_good, = axs[0].plot([], [], 'go')
ld_bad, = axs[0].plot([], [], 'ro', mfc='None')
lrs, = axs[0].plot([], [], 'bo-')

axs[1].set_ylim(0, 45)
axs[1].set_xlabel('Range Cell')
axs[1].set_ylabel('Monopole (A3) SNR (dB)')
ls_good, = axs[1].plot([], [], 'go')
ls_bad, = axs[1].plot([], [], 'ro', mfc='None')

# fig.legend((ld_good,ld_bad, lrs), ('good', 'badflagged', 'wtd averge'), 'upper left')

def subset_rsdata(d, c, bearing):
    # get data from qc'd and averaged (now data for RadialShorts) array
    xrow = numpy.where( (d[:,c['BEAR']]==bearing) & (d[:,c['VFLG']]==0) )[0]
    xcol = numpy.array([c['VELO'], c['SPRC'], c['BEAR']])
    a = d[numpy.ix_(xrow, xcol)]
    return a

def subset_data_good(d, c, bearing):
    # get GOOD data from RadialMetric array that is not badflagged 
    xrow = numpy.where( (d[:,c['BEAR']]==bearing) & (d[:,c['VFLG']]==0) )[0]
    xcol = numpy.array([c['VELO'], c['SPRC'], c['BEAR'], c['MA3S'], c['MSEL'], c['MSP1'], c['MDP1'], c['MDP2']])
    a = d[numpy.ix_(xrow, xcol)]
    return a

def subset_data_bad(d, c, bearing):
    # get BAD data from RadialMetric array that is badflagged 
    xrow = numpy.where( (d[:,c['BEAR']]==bearing) & (d[:,c['VFLG']]>0) )[0]
    xcol = numpy.array([c['VELO'], c['SPRC'], c['BEAR'], c['MA3S'], c['MSEL'], c['MSP1'], c['MDP1'], c['MDP2']])
    a = d[numpy.ix_(xrow, xcol)]
    return a

def init_plot(d, types_str, rsd, rstypes_str):
    """ 
    """
    c = get_columns(types_str)
    rsc = get_columns(rstypes_str)
    allranges = numpy.unique(d[:,c['SPRC']])
    allbearings = numpy.unique(d[:,c['BEAR']])
    bearing = allbearings[0]
    
    axs[0].set_xlim(allranges.min(), allranges.max())
    axs[1].set_xlim(allranges.min(), allranges.max())

    rs = subset_rsdata(rsd, rsc, bearing)
    gd = subset_data_good(d, c, bearing)
    bd = subset_data_bad(d, c, bearing)

    ld_good.set_xdata = gd[:,1]
    ld_good.set_ydata = gd[:,0]

    ld_bad.set_xdata = bd[:,1]
    ld_bad.set_ydata = bd[:,0]

    lsr.set_xdata = rs[:,1]
    lsr.set_ydata = rs[:,0]


def get_data(datadir, fn, patterntype, weight_parameter='MP'):
    """
    """
    # read in the data
    ifn = os.path.join(datadir, 'RadialMetric', patterntype, fn)
    d, types_str, header, footer = read_lluv_file(ifn)

    thresholds = [5.0, 50.0, 5.0]
    numfiles = 3
    numdegrees = 3
    weight_parameter = 'MP'

    ixfns = find_files_to_merge(ifn, numfiles, sample_interval=30)
    for xfn in ixfns:
        if xfn == ifn:
            continue
        d1, types_str1, _, _ = read_lluv_file(xfn)
        if len(d.shape) == len(d1.shape) == 2:
            if (d.shape[1] == d1.shape[1]) & (types_str == types_str1):
                # if same number and order of columns as d, then append the data d
                print '... ... merging: %s' % xfn
                d = numpy.vstack((d,d1))

    # do any threshold qc
    d = threshold_qc_all(d, types_str, thresholds)
   
    # do weighted averaging
    rsd, rstypes_str = weighted_velocities(d, types_str, numdegrees, weight_parameter)
    
    return d, types_str, rsd, rstypes_str


datadir = '/Users/codar/Documents/reprocessing_2015/Reprocess_HATY_70_35/'
# patterntype = 'MeasPattern' 
patterntype = 'IdealPattern' 
fn = 'RDLv_HATY_2013_11_05_0000.ruv'
d, types_str, rsd, rstypes_str = get_data(datadir, fn, patterntype, weight_parameter='MP')
    
