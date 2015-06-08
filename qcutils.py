#!/usr/bin/env python
#
# Last modified: Time-stamp: <2015-06-05 17:36:23 Sara>

"""Quality control (QC) functions for CODAR SeaSonde Radialmetric data

QC categories:
A. Threshold Tests -- badflag any values that fall below or above a single threshold
B. Range Tests -- badflag values that fall outside of a range
C. Weighted Averaging -- average several values with weights based on signal quality parameters 

QC Range Tests:
NONE TO DO YET

QC Threshold Tests:
1. DOA peak power (MSR1, MDR1, MDR2) < 5 dB default 
2. DOA 1/2 power width (3dB down) (MSW1, MDW1, MDW2) > 50 deg default
3. SNR on monopole (MA3S) < 5 dB default

Weighted Averaging:
1. Weighting based on Music Power (MSP1, MDP1, MDP2)
2. Weighting based on SNR on monopole (MA3S)
3. No weight function (None) 

"""
import sys
import os
import re
import fnmatch
import datetime

import numpy
numpy.set_printoptions(suppress=True)

# 
from codarutils import *

def _commonly_assigned_columns():
    """
    Commonly assigned CODAR RadialMetric columns
    """
    # for cut and pasting to other functions
    # help make some code more readable
    VFLG = c['VFLG']
    MSEL = c['MSEL']
    MSR1 = c['MSR1']
    MDR1 = c['MDR1']
    MDR2 = c['MDR2']
    MSW1 = c['MSW1'] 
    MDW1 = c['MDW1']
    MDW2 = c['MDW2']
    MA3S = c['MA3S']

def threshold_qc_doa_peak_power(d, types_str, threshold=5.0):
    """Bad Flag any DOA peak power (dB) less than threshold value (default 5.0 dB).

    Flags any direction of arrival (DOA) peak power (dB) that falls
    below the input threshold value (default 5.0 dB).  Depending on
    the value of MSEL (1, 2, or 3), MSR1, MDR1, or MDR2 columns are
    evaluated.  Returns modified matrix with VFLG column the only
    changed values.

    """
    c = get_columns(types_str)
    VFLG = c['VFLG'] # help make the test more readable
    MSEL = c['MSEL']
    MSR1 = c['MSR1']
    MDR1 = c['MDR1']
    MDR2 = c['MDR2']
    d1=numpy.copy(d)
    # 
    havenan = numpy.isnan(d[:,MSR1]) | numpy.isnan(d[:,MDR1]) | numpy.isnan(d[:,MDR2])
    bad = ((d[:,MSEL]==1) & (d[:,MSR1]<float(threshold))) | \
          ((d[:,MSEL]==2) & (d[:,MDR1]<float(threshold))) | \
          ((d[:,MSEL]==3) & (d[:,MDR2]<float(threshold))) | havenan
    d1[bad, VFLG] = d[bad,VFLG]+(1<<1)
    return d1

def threshold_qc_doa_half_power_width(d, types_str, threshold=50.0):
    """Bad Flag DOA 1/2 Power Width (degress) greater than threshold value (default 50.0 degrees).

    Flags any direction of arrival (DOA) 1/2 Power width (degress)
    that is wider than the input threshold value (default 50.0
    degrees).  Depending on the value of MSEL (1, 2, or 3), MSW1,
    MDW1, or MDW2 columns are evaluated.  Returns modified matrix with
    VFLG column the only changed values.

    """
    c = get_columns(types_str)
    VFLG = c['VFLG'] # help make the test more readable
    MSEL = c['MSEL']
    MSW1 = c['MSW1'] 
    MDW1 = c['MDW1']
    MDW2 = c['MDW2']
    d2=numpy.copy(d)
    # 
    havenan = numpy.isnan(d[:,MSW1]) | numpy.isnan(d[:,MDW1]) | numpy.isnan(d[:,MDW2])
    bad = ((d[:,MSEL]==1) & (d[:,MSW1]>float(threshold))) | \
          ((d[:,MSEL]==2) & (d[:,MDW1]>float(threshold))) | \
          ((d[:,MSEL]==3) & (d[:,MDW2]>float(threshold))) | havenan
    d2[bad, VFLG] = d[bad,VFLG]+(1<<2)
    return d2

def threshold_qc_monopole_snr(d, types_str, threshold=5.0):
    """Bad flag any SNR on monopole (dB)  less than threshold value (default 5.0 dB).

    Flags any signal-to-noise ratio (SNR) on monopole (dB) that falls
    below the input threshold value (default 5.0 dB) for all MSEL selections.

    """
    # Test 3 SNR on monopole (dB) for all selections
    # 
    c = get_columns(types_str)
    VFLG = c['VFLG'] # help make the test more readable
    MA3S = c['MA3S']
    d3=numpy.copy(d)
    # 
    bad = d[:,MA3S]<float(threshold)
    d3[bad, VFLG] = d[bad,VFLG]+(1<<3)
    return d3

def threshold_qc_all(d, types_str, thresholds=[5.0, 50.0, 5.0]):
    """Combine all three threshold tests

    Returns modified matrix with VFLG column only changed values.

    """
    # TO DO: Use a dict or list to pass thresholds for all tests ??
    # if dict what key to use
    # if list how know the correct order for tests
    # 
    dall = threshold_qc_doa_peak_power(d, types_str, thresholds[0])
    dall = threshold_qc_doa_half_power_width(dall, types_str, thresholds[1])
    dall = threshold_qc_monopole_snr(dall, types_str, thresholds[2])
    return dall


def weighted_velocities(d, types_str, bearing_spread=1.0, weight_parameter='MP'):
    """Calculates weighted average of radial velocities (VELO) at bearing and range.

    The weighted average of velocities found at given range and
    bearing based on weight_parameter.

    Paramters
    ---------
    d : ndarray
        The data from LLUV file(s). 
    types_str : string 
        The 'TalbleColumnTypes' string header of LLUV file(s) provide keys for each column.
    weight_parameter : string ('MP', 'SNR3', 'NONE'), optional 
        If 'MP' (default), uses MUSIC antenna peak power values for weighting function
           using MSEL to select one of (MSP1, MDP1, or MDP2).
        If 'SNR3', uses signal-to-noise ratio on monopole (MA3S).
        If 'NONE', just average with no weighting performed.
    bearing_spread: float, optional (default 1.0 degree)
       The number of degrees to look adjacent to either side of
       current bearing for increase spatial coverage (and potentially
       the sample size) for the average.
       For example, 
          If 0.0, velocities from window of 1 deg will be averaged.
          If 1.0, velocities from a window of 3 degrees will be avearged. This is the default.
          If 2.0, velocities from a window of 5 degrees

    Returns
    -------
    xd : ndarray
       The averaged values with range and bearing.
       An array with averaged values, range, bearing, 
    xtypes_str : string 
        The order and key-labels for each column of xd array

    """
    c = get_columns(types_str)
    offset = bearing_spread # 0 is default but can use 1 or 2 to get 3 or 5 deg spread
    # 
    ud = unique_rows(d[:,[c['SPRC'],c['BEAR'],c['VFLG']]].copy())
    # return only rows that have VFLG==0 (0 == good, >0 bad) so only get good data
    ud = ud[ud[:,2]==0]
    # sort this array based on rangecell (SPRC) and bearing (BEAR), remember last (col=0) is first to sort 
    idx = numpy.lexsort((ud[:,1], ud[:,0]))
    ud = ud[idx,:]
    # 
    # order of columns and labels for output data
    xtypes_str = 'SPRC BEAR VELO ESPC MAXV MINV EDVC ERSC'
    xc = get_columns(xtypes_str)
    #
    nrows, _ = ud.shape
    ncols = len(xc)
    xd = numpy.ones(shape=(nrows,ncols))*numpy.nan
    #
    for irow, cell in enumerate(ud):
        rngcell, bearing = cell[0:2]
        # numpy.where() returns a tuple for condition so use numpy.where()[0]
        # also VFLG must equal 0 (0 == good, >0 bad) so only get good data
        xrow = numpy.where((d[:,c['SPRC']]==rngcell) & \
                           (d[:,c['BEAR']]>=bearing-offset) & \
                           (d[:,c['BEAR']]<=bearing+offset) & \
                           (d[:,c['VFLG']]==0))[0]
        # If no row matches rngcell AND bearing, then no VELO data, skip to next bearing
        if xrow.size == 0: 
            continue
        
        xcol = numpy.array([c['VELO'], c['MSEL'], c['MSP1'], c['MDP1'], c['MDP2'], c['MA3S']])
        a = d[numpy.ix_(xrow, xcol)].copy()
        # if xrow.size == edvc:
        VELO = a[:,0] # all radial velocities found in cell
        SNR3 = a[:,5] # SNR on monopole for each velocity
        if weight_parameter.upper() == 'MP':
            # Create array to hold each Music Power (based on MSEL)
            MP = numpy.array(numpy.ones(VELO.shape)*numpy.nan) 
            # pluck the msel-based Music Power from MSP1, MDP1 or MPD2 column
            for msel in [1, 2, 3]:
                which = a[:,1]==msel
                MP[which,] = a[which, msel+1]
            # convert MP from db to voltage for weighting
            MP = numpy.power(10, MP/10.)
            wts = MP/MP.sum()
            velo = numpy.dot(VELO,wts)
        elif weight_parameter.upper() == 'SNR3' or weight_parameter.upper() == 'SNR':
            wts = SNR3/SNR3.sum()
            velo = numpy.dot(VELO,wts)
        elif weight_parameter.upper() == 'NONE':
            # do no weighting and just compute the mean of all velo's
            velo = VELO.mean()
        # data
        xd[irow,xc['SPRC']] = rngcell
        xd[irow,xc['BEAR']] = bearing
        xd[irow,xc['VELO']] = velo
        # other stat output
        xd[irow,xc['ESPC']] = VELO.std() # ESPC
        xd[irow,xc['MAXV']] = VELO.max() # MAXV
        xd[irow,xc['MINV']] = VELO.min() # MINV
        # (EDVC and ERSC are the same in this subroutine's context)
        xd[irow,xc['EDVC']] = VELO.size # EDVC Velocity Count 
        xd[irow,xc['ERSC']] = VELO.size # ERSC Spatial Count
                
    return xd, xtypes_str


def recursive_glob(treeroot, pattern):
    """ Glob-like search for filenames based on pattern but in all
    subdirectories to treeroot.

    Parameters
    ----------
    treeroot : string
       The top most directory path to begin search.
    pattern : string
       The pattern to match file in search.

    Return
    ------
    results : list of paths
       The results of search.

    >>> files = os.path.join(os.path.curdir, 'test', 'files')
    >>> recursive_glob(files, 'RDLx*.*')
    ['.\\test\\files\\codar_raw\\Radialshorts_HATY_2013_11_05\\RDLx_HATY_2013_11_05_0000.ruv',
    '.\\test\\files\\RadialShorts_mp_weight_angres1\\RDLx_HATY_2013_11_05_0000.ruv',
    '.\\test\\files\\RadialShorts_mp_weight_angres3\\RDLx_HATY_2013_11_05_0000.ruv',
    '.\\test\\files\\RadialShorts_no_weight_angres1\\RDLx_HATY_2013_11_05_0000.ruv',
    '.\\test\\files\\RadialShorts_snr_weight_angres1\\RDLx_HATY_2013_11_05_0000.ruv',
    '.\\test\\files\\RadialShorts_snr_weight_angres3\\RDLx_HATY_2013_11_05_0000.ruv']
    
    """
    # fnmatch gives you exactly the same patterns as glob, so this is
    # really an excellent replacement for glob.glob with very close
    # semantics.  Pasted from
    # <http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python>
    results = [] 
    for base, dirs, files in os.walk(treeroot): 
        goodfiles = fnmatch.filter(files, pattern) 
        results.extend(os.path.join(base, f) for f in goodfiles) 
    return results 

def _find_files():
    # one of the utilities for qc codar is to aggregate several sample times of data before qc'ing and doing
    # weighted average
    #
    # if interested in loading a targe site, date and time e.g. RDLv_HATY_2013_11_05_0000.ruv,
    # then we nned to first find it and the other
    # files (one before and one after) if they exist.
    # Need to know OutputTimeInterval (30 min in this case) to know what was previous file(s) and
    # next file(s) look like
    #
    # Most likely they will be in same directory but they might not
    # 
    # test/files/codar_raw/Radialmetric_HATY_2013_11_04/RDLv_HATY_2013_11_04_2330.ruv 
    # test/files/codar_raw/Radialmetric_HATY_2013_11_05/RDLv_HATY_2013_11_05_0000.ruv 
    # test/files/codar_raw/Radialmetric_HATY_2013_11_05/RDLv_HATY_2013_11_05_0030.ruv
    files = recursive_glob(os.path.join(files, 'codar_raw'), 'RDLv*.ruv')
    for f in files:
        fn = os.path.split(f)[-1]
        fndt = filt_datetime(fn)

    target
    dt_start = dts[0]-datetime.timedelta(minutes=30)
    dt_end = dts[0]+datetime.timedelta(minutes=30)

def filt_datetime(input_string):
    """ Attempts to filter date and time from input string.
    
    Following the template, YYYY(-)MM(-)DD(-)(hh(:)(mm(:)(ss)))
    find the most precise, reasonable string match and
    return its datetime object.

    Typical matches include, YYYYMMDD-hhmmss, YYYY-MM-DD-hh:mm:ss

    Requires date with all three (year, month, day) in decreasing
    order as integers. Time is optional.
    
    """
    # typical codar time stamp format
    pattern = r"""
    # YYYY(-)MM(-)DD(-)(hh(:)(mm(:)(ss)))
    (\d{4})           # 4-digit YEAR 
    \D?               # optional 1 character non-digit separator (e.g. ' ' or '-')
    (\d{2})           # 2-digit MONTH 
    \D?               # optional 1 character non-digit separator
    (\d{2})           # 2-digit DAY 
    \D?               # optional 1 character non-digit separator (e.g. ' ' or 'T')
    (\d{2})?          # optional 2-digit HOUR 
    \D?               # optional 1 character non-digit separator (e.g. ' ' or ':')
    (\d{2})?          # optional 2-digit MINUTE 
    \D?               # optional 1 character non-digit separator (e.g. ' ' or ':')
    (\d{2})?          # optional 2-digit SECOND
    """
    p = re.compile(pattern, re.VERBOSE)
    # input_string = 'RDLv_HATY_2013_11_05_000000.ruv'
    # input_string = 'RDLv_HATY_2013_11_05_0000.ruv'
    # input_string = 'RDLv_HATY_2013_11_05_00.ruv'
    # input_string = 'RDLv_HATY_2013_11_05.ruv'
    # input_string = 'RDLv_HATY_13_11_05.ruv'
    m = p.search(input_string) 
    # m.groups() # should be ('2013', '11', '05', '00', '00', None) for 'RDLv_HATY_2013_11_05_0000.ruv'
    if m:
        values = [int(yi) for yi in m.groups() if yi is not None] # [2013, 11, 5, 0, 0]
        # datetime.datetime(*v) requires mininum of year, month, day
        dt = datetime.datetime(*values) # datetime.datetime(2013, 11, 5, 0, 0)
    else:
        dt = None
    return dt

# for testing
if __name__ == '__main__':
    # 
    # datadir = sys.argv[1]
    # patterntype = sys.argv[2]
    # patterntype = 'MeasPattern' 
    # patterntype = 'IdealPattern'
    
    # read in the data
    ifn = os.path.join('.', 'test', 'files', 'codar_raw', \
                   'Radialmetric_HATY_2013_11_05', \
                   'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    # thresholding
    # dall = threshold_qc_all(d, types_str, thresholds=[5.0, 50.0, 5.0])
    # weighting
    xd, xtypes_str = weighted_velocities(d, types_str, 0.0, 'SNR')
    
    # create radialshort data, first generate array then fill it
    rsd, rsdtypes_str = generate_radialshort_array(d, types_str)
    rsd = fill_radialshort_array(rsd, rsdtypes_str, xd, xtypes_str)[0]
    # modify header from radialmetric, based on new radialshort data
    rsdheader = generate_radialshort_header(rsd, rsdtypes_str, header)
    # not modifying the footer at this time
    rsdfooter = footer
    # output the radialshort data specified location
    ofn = os.path.join('.', 'test', 'files', 'test_output.txt')
    write_output(ofn, rsdheader, rsd, rsdfooter)
    
    #
    rsc = get_columns(rsdtypes_str)
    xc = get_columns(xtypes_str)
