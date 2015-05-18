"""
QC functions for radialmetric data


"""
def _asssign_columns():
    # for cut and pasting to other functions
    # help make the test more readable
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

    Returns modified matrix with VFLG column the only changed values.

    """
    # 
    dall = threshold_qc_doa_peak_power(d, types_str, thresholds[0])
    dall = threshold_qc_doa_half_power_width(dall, types_str, thresholds[1])
    dall = threshold_qc_monopole_snr(dall, types_str, thresholds[2])
    return dall


def weighted_velocities(d, type_str, xd, xtype_str, bearing_spread=0, weight_parameter='MP'):
    # 
    c = get_columns(type_str)
    xc = get_columns(xtype_str)
    d1=numpy.copy(d)
    offset = bearing_spread # 0 is default but can use 1 or 2 to get 3 or 5 deg spread
    # 
    (nrows, ncols) = d.shape
    for irow in range(nrows):
        rngcell = d[irow,c['SPRC']]
        bearing = d[irow,c['BEAR']]
        vflg = d[irow,c['VFLG']]
        edvc = d[irow,c['EDVC']]
        head = d[irow,c['HEAD']]
        # print "%d %d %d %d %d %5.1f" % (irow, rngcell, bearing, vflg, edvc, head)
        #
        # if not flagged in radial shorts, recompute VELO, and (VELU and VELV) based on weighting from radial metric
        if not vflg:
            # find rows in radial metric data at same range and bearing, where VLFG also is 0
            # numpy.where() return a tuple for condition so use numpy.where()[0]
            # xrow = numpy.where((xd[:, xc['SPRC']]==rngcell) & (xd[:,xc['BEAR']]==bearing) & (xd[:,xc['VFLG']]==0))[0]
            xrow = numpy.where((xd[:, xc['SPRC']]==rngcell) & (xd[:,xc['BEAR']]>=bearing-offset) & (xd[:,xc['BEAR']]<=bearing+offset) & (xd[:,xc['VFLG']]==0))[0] 
            xcol = numpy.array([xc['VELO'], xc['MSEL'], xc['MSP1'], xc['MDP1'], xc['MDP2'], xc['MA3S']])
            a = xd[numpy.ix_(xrow, xcol)].copy()
            # if xrow.size == edvc:
            VELO = a[:,0] # all radial velocities found in cell
            SNR3 = a[:,5] # SNR on monopole for each velocity

            if weight_parameter == 'MP':
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
            elif weight_parameter == 'SNR3':
                wts = SNR3/SNR3.sum()
                velo = numpy.dot(VELO,wts)
            elif weight_parameter.upper() == 'NONE':
                # do no weighting and just compute the mean of all velo's
                velo = VELO.mean()
 
            (velu, velv) = compass2uv(velo, head)
            # print '... %10.3f %5.1f %10.3f %10.3f' % (velo, head, velu, velv)
            # VELO.mean()
            # replace VELO, VELU, VELV in radial short data
            d1[irow,c['VELO']]=velo
            d1[irow,c['VELU']]=velu
            d1[irow,c['VELV']]=velv
                
    return d1
