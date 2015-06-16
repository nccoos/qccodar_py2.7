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
import matplotlib.pylab as pylab

# initialize the main parameters that can be changed on the plot
params = {'thresholds' : [5.0, 50.0, 5.0],
          'numfiles' : 3,
          'numdegrees' :  3,
          'weight_parameter' : 'SNR',
          'bearing' : 0,
}

pylab.rcParams['figure.figsize']= (11.0, 8.0)
fig, axs = plt.subplots(3,1)

axs[0].axhline(y=0, linewidth=1, color='k')
axs[0].set_ylim(-150, 150)
axs[0].set_ylabel('Radial Velocity (cm/s)')
ld_good, = axs[0].plot([], [], 'go', mec='g')
ld_bad, = axs[0].plot([], [], 'ro', mec='r', mfc='None')
lrs, = axs[0].plot([], [], 'bo-')

axs[1].set_ylim(0, 45)
axs[1].set_xlabel('Range Cell')
axs[1].set_ylabel('Monopole (A3) SNR (dB)')
ls_good, = axs[1].plot([], [], 'go', mec='g')
ls_bad, = axs[1].plot([], [], 'ro', mec='r', mfc='None')

# change position of last plot to stay on left margin but make it square
bb2 = axs[2].get_position()
bb2.bounds = (bb2.bounds[0], bb2.bounds[1], bb2.height, bb2.height)
axs[2].set_position(bb2)

axs[2].set_ylim(-1,1)
axs[2].set_xlim(-1,1)
axs[2].axhline(y=0, linewidth=1, color='k')
axs[2].axvline(x=0, linewidth=1, color='k')
axs[2].set_aspect('equal')
axs[2].set_xticklabels('')
axs[2].set_yticklabels('')

lbear, = axs[2].plot([0,compass2uv(1,45)[0]], [0,compass2uv(1,45)[1]], 'b-')

# Widget functions
def sbear_change(val):
    global params
    params['bearing'] = numpy.round(val)
    plot_data(d, types_str, rsd, rstypes_str)
    fig.canvas.draw()

def wtdavg_change(label):
    global params
    global rsd
    params['weight_parameter'] = label

    if params['weight_parameter']=='MP':
        axs[1].set_ylim(-125, -75)
        axs[1].set_ylabel('MUSIC Power (dB)')
    elif params['weight_parameter']=='SNR':
        axs[1].set_ylim(0, 45)
        axs[1].set_ylabel('Monopole (A3) SNR (dB)')
    elif params['weight_parameter']=='NONE':
        axs[1].set_ylim(0, 1)
        axs[1].set_ylabel('No Weighting Param')
 
    rsd, rstypes_str = weighted_velocities(d, types_str, params['numdegrees'], params['weight_parameter'])
    plot_data(d, types_str, rsd, rstypes_str)
    fig.canvas.draw()

# Widgets
axbear = plt.axes([0.1, 0.05, 0.8, 0.03]) 
sbear = matplotlib.widgets.Slider(axbear, 'Bearing', 0, 359, valinit=0, valfmt=u'%d')
sbear.on_changed(sbear_change)

axradio = plt.axes([0.4, 0.1, 0.15, 0.15], aspect='equal')
rwtdavg = matplotlib.widgets.RadioButtons(axradio, ('MP', 'SNR', 'NONE'), active=1)
rwtdavg.on_clicked(wtdavg_change)

# fig.legend((ld_good,ld_bad, lrs), ('good', 'badflagged', 'wtd averge'), 'upper left')

def subset_rsdata(d, c, bearing):
    # get data from qc'd and averaged (now data for RadialShorts) array
    xrow = numpy.where( (d[:,c['BEAR']]==bearing) )[0]
    xcol = numpy.array([c['VELO'], c['SPRC'], c['BEAR']])
    a = d[numpy.ix_(xrow, xcol)]
    return a

def subset_data_good(d, c, bearing):
    # get GOOD data from RadialMetric array that is not badflagged 
    xrow = numpy.where( (d[:,c['BEAR']]==bearing) & (d[:,c['VFLG']]==0) )[0]
    xcol = numpy.array([c['VELO'], c['SPRC'], c['BEAR'], c['MA3S'], c['MSEL'], c['MSP1'], c['MDP1'], c['MDP2']])
    a = d[numpy.ix_(xrow, xcol)]
    # Create array to hold each Music Power (based on MSEL)
    MP = numpy.array(numpy.ones(a[:,0].shape)*numpy.nan) 
    # pluck the msel-based Music Power from MSP1, MDP1 or MPD2 column
    for msel in [1, 2, 3]:
        which = a[:,4]==msel
        MP[which,] = a[which, msel+4]
    # append as last column
    a = numpy.hstack((a,MP.reshape(MP.size,1)))
    return a

def subset_data_bad(d, c, bearing):
    # get BAD data from RadialMetric array that is badflagged 
    xrow = numpy.where( (d[:,c['BEAR']]==bearing) & (d[:,c['VFLG']]>0) )[0]
    xcol = numpy.array([c['VELO'], c['SPRC'], c['BEAR'], c['MA3S'], c['MSEL'], c['MSP1'], c['MDP1'], c['MDP2']])
    a = d[numpy.ix_(xrow, xcol)]
    # Create array to hold each Music Power (based on MSEL)
    MP = numpy.array(numpy.ones(a[:,0].shape)*numpy.nan) 
    # pluck the msel-based Music Power from MSP1, MDP1 or MPD2 column
    for msel in [1, 2, 3]:
        which = a[:,4]==msel
        MP[which,] = a[which, msel+4]
    # append as last column
    a = numpy.hstack((a,MP.reshape(MP.size,1)))
    return a

def az2deg(az):
    x,y = compass2uv(1,az)
    # arctan2 does this relative to (x0,y0)=(1,0)
    return numpy.arctan2(y, x)*180./numpy.pi

def init_plot(d, types_str, rsd, rstypes_str):
    """ Set plot and slider limits
    """
    global params

    c = get_columns(types_str)
    rsc = get_columns(rstypes_str)
    xrow = numpy.where( (d[:,c['VFLG']]==0) )[0]
    allranges = numpy.unique(d[:,c['SPRC']] )
    allbearings = numpy.unique(d[xrow,c['BEAR']])
    params['bearing'] = allbearings[0]
    
    axs[0].set_xlim(0, allranges.max()+2)
    axs[1].set_xlim(0, allranges.max()+2)

    # add radar patch limits
    thetas = numpy.array([az2deg(b) for b in allbearings])
    axs[2].add_patch(matplotlib.patches.Wedge((0,0), 1, thetas.min(), thetas.max(), \
                                              zorder=-1, ec='None', fc=(.9,.9,.9)))

    # update the slider with initialized bearings
    sbear.valinit = params['bearing']
    sbear.valmin = allbearings[0]
    sbear.valmax = allbearings[-1]

    # put the first bearing data into the plots
    plot_data(d, types_str, rsd, rstypes_str)


def plot_data(d, types_str, rsd, rstypes_str):
    # print bearing

    # update new bearing line
    lbear.set_xdata([0, compass2uv(1,params['bearing'])[0]])
    lbear.set_ydata([0, compass2uv(1,params['bearing'])[1]])

    c = get_columns(types_str)
    rsc = get_columns(rstypes_str)

    # update plots with new bearing
    rs = subset_rsdata(rsd, rsc, params['bearing'])
    gd = subset_data_good(d, c, params['bearing'])
    bd = subset_data_bad(d, c, params['bearing'])

    velo = 0
    rnge = 1
    mp = 8
    snr = 3
    if gd.size>0:
        ld_good.set_xdata(gd[:,rnge])
        ld_good.set_ydata(gd[:,velo])
        if params['weight_parameter']=='MP':
            ls_good.set_xdata(gd[:,rnge])
            ls_good.set_ydata(gd[:,mp])
        elif params['weight_parameter']=='SNR':
            ls_good.set_xdata(gd[:,rnge])
            ls_good.set_ydata(gd[:,snr])
        elif params['weight_parameter']=='NONE':
            ls_good.set_xdata([])
            ls_good.set_ydata([])

    if bd.size>0:
        ld_bad.set_xdata(bd[:,rnge])
        ld_bad.set_ydata(bd[:,velo])
        if params['weight_parameter']=='MP':
            ls_bad.set_xdata(bd[:,rnge])
            ls_bad.set_ydata(bd[:,mp])
        elif params['weight_parameter']=='SNR':
            ls_bad.set_xdata(bd[:,rnge])
            ls_bad.set_ydata(bd[:,snr])
        elif params['weight_parameter']=='NONE':
            ls_bad.set_xdata([])
            ls_bad.set_ydata([])

    if rs.size>0:
        lrs.set_xdata(rs[:,1])
        lrs.set_ydata(rs[:,0])


def get_data(datadir, fn, patterntype):
    """
    """
    # read in the data
    ifn = os.path.join(datadir, 'RadialMetric', patterntype, fn)
    d, types_str, header, footer = read_lluv_file(ifn)

    ixfns = find_files_to_merge(ifn, params['numfiles'], sample_interval=30)
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
    d = threshold_qc_all(d, types_str, params['thresholds'])
   
    # do weighted averaging
    rsd, rstypes_str = weighted_velocities(d, types_str, params['numdegrees'], params['weight_parameter'])
    
    return d, types_str, rsd, rstypes_str


datadir = os.path.join('.', 'test', 'files', 'codar_raw')
# datadir = '/Users/codar/Documents/reprocessing_2015/Reprocess_HATY_70_35/'
# patterntype = 'MeasPattern' 
patterntype = 'IdealPattern' 
fn = 'RDLv_HATY_2013_11_05_0000.ruv'
d, types_str, rsd, rstypes_str = get_data(datadir, fn, patterntype)
init_plot(d, types_str, rsd, rstypes_str)
