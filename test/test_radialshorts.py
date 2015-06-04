#!/usr/bin/env python
#
# Last modified: Time-stamp: <2015-06-01 15:46:36 haines>
"""
Test functions for generating radialshort data output.

"""
import os
import numpy
numpy.set_printoptions(suppress=True)
from qcutils import *

files = os.path.join(os.path.curdir, 'test', 'files')

def test_generate_radialshort_array():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    ####################
    # TESTING
    rsd, rsdtypes_str = generate_radialshort_array(d, types_str) 
    ####################
    rsc = get_columns(rsdtypes_str)

    ifn2 = os.path.join(files, 'RadialShorts_no_weight_angres1', 'RDLx_HATY_2013_11_05_0000.ruv')
    td, ttypes_str, theader, tfooter = read_lluv_file(ifn2)
    tc = get_columns(ttypes_str)

    trngbear = td[:, [tc['SPRC'], tc['BEAR']]]
    rsrngbear = rsd[:, [rsc['SPRC'], rsc['BEAR']]]
    trows, rsrows = cell_intersect(trngbear, rsrngbear)
    # trows, rsrows = cell_intersect(td[:, [tc['SPRC'], tc['BEAR']]], rsd[:, [rsc['SPRC'], rsc['BEAR']]])
    # 
    subtd = td[trows, :]
    subrsd = rsd[rsrows, :]
    # these are the columns generated from d that we want to compare with CODAR radialshorts
    tcol  = numpy.array([tc['LOND'], tc['LATD'], tc['VFLG'], tc['RNGE'], tc['BEAR'], tc['SPRC']])
    rscol = numpy.array([rsc['LOND'], rsc['LATD'], rsc['VFLG'], rsc['RNGE'], rsc['BEAR'], rsc['SPRC']])
    # subrsd is close to subtd  within 1/1000 th, since test data was output by CODAR
    assert numpy.isclose(subrsd[:,rscol], subtd[:,tcol], rtol=1e-05, atol=1e-03, equal_nan=True).all()

def test_fill_radialshort_array():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    rsd, rsdtypes_str = generate_radialshort_array(d, types_str)
    
    # no qc thresholds and this no weighting and 1 deg angres (bearing_offset=0.0) is same as CODAR processing
    xd, xtypes_str = weighted_velocities(d, types_str, 0.0, 'NONE')
    ####################
    # TESTING
    rsd = fill_radialshort_array(rsd, rsdtypes_str, xd, xtypes_str)[0] # <----------TESTING
    ####################
    rsc = get_columns(rsdtypes_str)

    ifn2 = os.path.join(files, 'RadialShorts_no_weight_angres1', 'RDLx_HATY_2013_11_05_0000.ruv')
    td, ttypes_str, theader, tfooter = read_lluv_file(ifn2)
    tc = get_columns(ttypes_str)

    trngbear = td[:, [tc['SPRC'], tc['BEAR']]]
    rsrngbear = rsd[:, [rsc['SPRC'], rsc['BEAR']]]
    trows, rsrows = cell_intersect(trngbear, rsrngbear)
    # trows, rsrows = cell_intersect(td[:, [tc['SPRC'], tc['BEAR']]], rsd[:, [rsc['SPRC'], rsc['BEAR']]])
    # 
    subtd = td[trows, :]
    subrsd = rsd[rsrows, :]

    # 1st make sure generated columsn okay
    # these are the columns generated from d that we want to compare with CODAR radialshorts
    tcol  = numpy.array([tc['LOND'], tc['LATD'], tc['VFLG'], tc['RNGE'], tc['BEAR'], tc['SPRC']])
    rscol = numpy.array([rsc['LOND'], rsc['LATD'], rsc['VFLG'], rsc['RNGE'], rsc['BEAR'], rsc['SPRC']])
    # subrsd is close to subtd  within 1/1000 th, since test data was output by CODAR
    assert numpy.isclose(subrsd[:,rscol], subtd[:,tcol], rtol=1e-05, atol=1e-03, equal_nan=True).all(), \
        'think that generate_radialshort_array failed this time'
    
    # 2nd then test the columns that are filled 
    # these are the columns filled from xd VELO and that we want to compare with CODAR radialshorts
    tcol  = numpy.array([tc['VELO']])
    rscol = numpy.array([rsc['VELO']])
    assert numpy.isclose(subrsd[:,rscol], subtd[:,tcol], rtol=1e-05, atol=1e-03, equal_nan=True).all(), \
        'something wrong with VELOs'

    tcol  = numpy.array([tc['XDST'], tc['YDST']])
    rscol = numpy.array([rsc['XDST'], rsc['YDST']])
    assert numpy.isclose(subrsd[:,rscol], subtd[:,tcol], rtol=1e-05, atol=1e-03, equal_nan=True).all(), \
        'something wrong with XDST, YDST'

    tcol  = numpy.array([tc['VELU'], tc['VELV']])
    rscol = numpy.array([rsc['VELU'], rsc['VELV']])
    assert numpy.isclose(subrsd[:,rscol], subtd[:,tcol], rtol=1e-05, atol=1e-03, equal_nan=True).all(), \
        'something wrong with VELU, VELV'

def test_generate_radialshort_header():
    pass

def _scratch():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    xd, xtypes_str = weighted_velocities(d, types_str, bearing_spread=0.0, weight_parameter='NONE')
    xc = get_columns(xtypes_str)

    ifn2 = os.path.join(files, 'RadialShorts_no_weight_angres1', 'RDLx_HATY_2013_11_05_0000.ruv')
    td, ttypes_str, theader, tfooter = read_lluv_file(ifn2)
    tc = get_columns(ttypes_str)

    trngbear = td[:, [tc['SPRC'], tc['BEAR']]]
    xrngbear = xd[:, [xc['SPRC'], xc['BEAR']]]
    trows, xrows = cell_intersect(trngbear, xrngbear)
    # trows, xrows = cell_intersect(td[:, [tc['SPRC'], tc['BEAR']]], xd[:, [xc['SPRC'], xc['BEAR']]])

    subtd = td[trows, tc['VELO']]
    subxd = xd[xrows, xc['VELO']]
    # subxd VELO is close to subtd VELO within 1/1000 th, since test data was output by CODAR
    assert numpy.isclose(subxd, subtd, rtol=1e-05, atol=1e-03, equal_nan=True).all()
