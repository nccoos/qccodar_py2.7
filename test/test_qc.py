#!/usr/bin/env python
#
# Last modified: Time-stamp: <2015-06-01 15:46:36 haines>
"""
Tests for qc thresholds and weighted averaging.

"""
import os
import numpy
numpy.set_printoptions(suppress=True)
from qcutils import *

files = os.path.join(os.path.curdir, 'test', 'files')

def test_read_test0():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    #
    d0 = d
    #
    ifn2 = os.path.join(files, 'Radialmetric_test0', 'RDLv_HATY_2013_11_05_0000.ruv')
    td, ttypes_str, theader, tfooter = read_lluv_file(ifn2)
    assert numpy.isclose(d0, td, equal_nan=True).all(), 'should be equal, including where NaN'

def test_threshold_qc_doa_peak_power():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    #
    d1 = threshold_qc_doa_peak_power(d, types_str)
    #
    ifn2 = os.path.join(files, 'Radialmetric_test1', 'RDLv_HATY_2013_11_05_0000.ruv')
    td, ttypes_str, theader, tfooter = read_lluv_file(ifn2)
    #
    assert numpy.isclose(d1, td, equal_nan=True).all(), 'should be equal, including where NaN'

def test_threshold_qc_doa_half_power_width():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    #
    d2 = threshold_qc_doa_half_power_width(d, types_str)
    #
    ifn2 = os.path.join(files, 'Radialmetric_test2', 'RDLv_HATY_2013_11_05_0000.ruv')
    td, ttypes_str, theader, tfooter = read_lluv_file(ifn2)
    #
    assert numpy.isclose(d2, td, equal_nan=True).all(), 'should be equal, including where NaN'

def test_threshold_qc_monopole_snr():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    #
    d3 = threshold_qc_monopole_snr(d, types_str)
    #
    ifn2 = os.path.join(files, 'Radialmetric_test3', 'RDLv_HATY_2013_11_05_0000.ruv')
    td, ttypes_str, theader, tfooter = read_lluv_file(ifn2)
    #
    assert numpy.isclose(d3, td, equal_nan=True).all(), 'should be equal, including where NaN'

def test_threshold_qc_all():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    #
    dall = threshold_qc_all(d, types_str)
    #
    ifn2 = os.path.join(files, 'Radialmetric_testall', 'RDLv_HATY_2013_11_05_0000.ruv')
    td, ttypes_str, theader, tfooter = read_lluv_file(ifn2)
    #
    assert numpy.isclose(dall, td, equal_nan=True).all(), 'should be equal, including where NaN'



def _scratch():
    ofn = os.path.join(files, 'test1_output.txt')
    write_output(ofn, header, d1, footer)
    #
    idx = numpy.where(d1 != td)    
    assert numpy.isnan(d[idx]).all()
    assert numpy.isnan(d2[idx]).all()
    for i,j in numpy.array(idx).T:
        # if not numpy.isnan(d1[i,j]):
            print "(%4d, %4d) %5g %5g" % (i,j, d1[i,j], td[i,j])    
