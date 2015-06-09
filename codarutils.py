#!/usr/bin/env python
# 
# Last modified: Time-stamp: <2015-06-05 17:05:49 Sara>
""" CODAR Utilities 

"""
import sys
import os
import re
import fnmatch
import datetime

import numpy
numpy.set_printoptions(suppress=True)
from StringIO import StringIO

def load_data(inFile):
    lines=None
    if os.path.exists(inFile):
        f = open(inFile, 'r')
        lines = f.readlines()
        f.close()
        if len(lines)<=0:
            print 'Empty file: '+ inFile           
    else:
        print 'File does not exist: '+ inFile
    return lines

def _write_empty_output(ofn, header, footer):
    """deprecating -- write header and footer only, since no radial data"""
    f = open(ofn, 'w')
    f.write(header)
    f.write(footer)
    f.close()

def write_output(ofn, header, d, footer):
    """Write header, radialmetric data, and footer. """
    f = open(ofn, 'w')
    if header[-1] == '\n':
        f.write(header)
    else:
        f.write(header+'\n')
    # if there is any data, save to the file)
    if d.size > 0:
        numpy.savetxt(f, d, fmt='%g')
    f.write(footer)
    f.close()

def read_lluv_file(ifn):
    """Reads header, CSV table, and tail of LLUV files.  

    Extracts LLUV data into numpy array for further processing. If
    there is radial data in the file, these lines are found in the
    middle and each line has no '%'. All header and footer lines start
    with '%'. The middle is bracketed by header and footer lines. This
    routine searches for a header, middle, and footer.  If no midde is
    found, it is assumed that no radial data exists for the site and
    time.  If no middle, all the comment lines fall in the header and
    the footer is empty.

    Parameter
    ---------
    ifn : string
       The input filename and path.

    Returns
    -------
    d : ndarray
       The radial data bound by header and footer.  If there is no data, 
       d is an empty array (d.size==0), then no radial table was found.
    types_str : string 
       The order and label of columns in d array.  If there is no data,
       types_str is an empty string ('').
    header : string 
       All the '%' commented lines preceding '%TableStart:'
    footer : string
       All the '%' commented lines after the data table, starting with '%TableEnd:'

    """
    lines = load_data(ifn)
    m=re.match(r'(?P<header>(%.*\n)*)(?P<middle>([\d\s-].*\n)*)(?P<tail>(%.*\n)*)', \
               ''.join(lines))
    header  = m.group('header')
    footer = m.group('tail')

    # did not find a middle, so all comments are in header, and footer is empty
    if len(footer)<=0:
        m = re.findall(r'^(%.*):\s*(.*)$', header, re.MULTILINE)
        for k,v in m:
            if k == '%TableColumnTypes':
                types_str = v
                break            
      
        print 'No Radial Data in '+ ifn
        return numpy.array([]), types_str, header, footer

    # read header that match '%(k): (v)\n' pairs on each line
    m = re.findall(r'^(%.*):\s*(.*)$', header, re.MULTILINE)
    for k,v in m:
        ### print k+', '+v
        if k == '%TimeStamp':
            #sample_dt = scanf_datetime(v, fmt='%Y %m %d %H %M %S')
            pass
        elif k == '%TableType':
            ftype = v
        elif k == '%TableColumns':
            ncol = int(v)
        elif k == '%TableRows':
            nrow = int(v)
        elif k == '%TableColumnTypes':
            types_str = v
        elif k == '%TableStart':
            break

    # use file object from lines to extract 
    s = StringIO(''.join(lines))
    s.seek(0) # ensures start posn of file-like string s
    d = numpy.loadtxt(s, comments='%')
    # lat, lon, u, v = numpy.loadtxt(s, usecols=(0,1,2,3), comments='%', unpack=True)
    return d, types_str, header, footer

def get_columns(types_str):
    # use dict to store column label and it's column number
    #c = col.defaultdict(int)
    c = {}
    column_labels = types_str.strip().split(' ')
    m = re.findall(r'\w{4}', types_str)
    for label in column_labels:
        c[label]=m.index(label) # c['VFLG']=4
    return c

def generate_radialshort_array(d, types_str, table_type='LLUV RDL7'):
    """Generates radialshort (rsd) data array.

    This function generates radialshort data (rsd) array by merging
    unique rows of original radialmetric data (d).  It is intended to be
    the output CSV data table (middle) of CODAR LLUV format with header,
    middle, and footer.
    
    Parameters
    ----------
    d : ndarray
       The original radialmetric data for extracting lat, lon, etc.
    types_str : string 
        The order and key-labels for each column of xd array.
    table_type : string

    Returns
    -------
    rsd : ndarray
       The generated radialshort (rsd) data from merging unique rows
       of rangecell and bearing.
    rsdtypes_str : string 
        The order and key-labels for each column of rsd array.

    """
    if table_type == 'LLUV RDL7':
        rsdtypes_str = 'LOND LATD VELU VELV VFLG ESPC MAXV MINV EDVC ERSC XDST YDST RNGE BEAR VELO HEAD SPRC'
    else:
        print 'generate_radial_array() : Unrecognized table_type "%s"' % (table_type,)
        return numpy.array([]), ''

    if d.size == 0:
        return numpy.array([]), rsdtypes_str

    # 
    # get unique rows of rangecells and bearings based on input radialmetric data
    # but also collect columns of info so we don't have to reproduce in next steps
    c = get_columns(types_str)
    dcol = numpy.array([c['LOND'], c['LATD'], c['VFLG'], c['RNGE'], c['BEAR'], c['SPRC']])
    ud = unique_rows(d[:,dcol].copy())
    # return only rows that have VFLG==0 (0 == good, >0 bad) so only get good data
    ud = ud[ud[:,2]==0]
    # sort this array based on rangecell (SPRC) and bearing (BEAR)
    idx = numpy.lexsort((ud[:,4], ud[:,5]))
    ud = ud[idx,:]
    # 
    # order of columns and labels for output data
    rsc = get_columns(rsdtypes_str)
    nrows,_ = ud.shape
    ncols = len(rsc)
    # Initialize new array for radial shorts
    rsd = numpy.ones(shape=(nrows,ncols))*numpy.nan
    # Distribute data in to new array 
    rscol = numpy.array([rsc['LOND'], rsc['LATD'], rsc['VFLG'], rsc['RNGE'], rsc['BEAR'], rsc['SPRC']])
    rsd[:,rscol] = ud
    return rsd, rsdtypes_str
    
def fill_radialshort_array(rsd, rsdtypes_str, xd, xtypes_str):
    """Fill in radialshort (rsd) data. 

    This function fills radialshort data (rsd) by merging rows of
    averaged velocity data (xd) where rangecell and bearing match.

    Parameters
    ----------
    rsd : ndarray
       The generated radialshort (rsd) data from generate_radialshort_array().
    rsdtypes_str : string 
        The order and key-labels for each column of rsd array.
    xd : ndarray
       The array with averaged values, range, bearing, and other stats used to fill in rsd array.
    xtypes_str : string 
        The order and key-labels for each column of xd array.

    Returns
    -------
    rsd : ndarray
       The modified radialshort (rsd) data.
    rsdtypes_str : string 
        The order and key-labels for each column of rsd array.

    """

    rsc = get_columns(rsdtypes_str)
    rscells = rsd[:,[rsc['SPRC'], rsc['BEAR']]]
    rscol = numpy.array([rsc['VELO'], rsc['ESPC'], rsc['MAXV'], rsc['MINV'], rsc['EDVC'], rsc['ERSC']])

    xc = get_columns(xtypes_str)
    xcells = xd[:,[xc['SPRC'],xc['BEAR']]]
    xcol = numpy.array([xc['VELO'], xc['ESPC'], xc['MAXV'], xc['MINV'], xc['EDVC'], xc['ERSC']])

    # check range and bearing columns are the same between rsd and xd
    assert rscells.shape == xcells.shape, "rscells.shape(%d,%d)" % rscells.shape 
    assert (rscells == xcells).all()
    # deal xd data into rsd by columns
    rsd[:,rscol] = xd[:,xcol]

    # if not, then match row-by-row, and deal each row as matches are made
    # for rngcell, bearing in rsd[:10, [rsc['SPRC'],rsc['BEAR']]]:
    #     print "rangecell: %d, bearing: %d" % (rngcell, bearing)

    # create HEAD column based on BEAR+180
    bear = rsd[:,rsc['BEAR']]
    head = numpy.mod(bear+180., 360.)
    velo = rsd[:,rsc['VELO']]
    # compute velocity components
    (velu, velv) = compass2uv(velo, head)
    # replace VELU, VELV, HEAD in radial short data
    rsd[:,rsc['VELU']]=velu
    rsd[:,rsc['VELV']]=velv
    rsd[:,rsc['HEAD']]=head
    #
    rnge = rsd[:,rsc['RNGE']]
    (xdist, ydist) = compass2uv(rnge,bear)
    # relace XDIST, YDIST in radial short data
    rsd[:,rsc['XDST']]=xdist
    rsd[:,rsc['YDST']]=ydist
    
    return rsd, rsdtypes_str

def generate_radialshort_header(rsd, rsdtypes_str, header):
    """ Fill radialshort header details from radialmetric header

    Replaces lines from input radialmetric header that start with '%
    Table*' and '%% ' to match radialshort format.  Information gleaned
    from rsd and rdstypes_str are used in header metadata.

    Note: this could be generalized for any LLUV file since passing in info
    create lines that describe the main data. LLUV RDL7 is specified in
    generate_radialshort_array() as is the rsdtypes_str defining columns.

    Parameter:
    ----------
    rsd : ndarray
       The radialshort (rsd) data.
    rsdtypes_str : string 
        The order and key-labels for each column of rsd array.
    header : string
       The radialmetric header

    Returns:
    --------
    rsdheader : string
       The radialshort header
    """

    # keep everything up until TableType
    rsdheader = re.split(r'(\n%TableType)', header)[0]

    ncols_from_string = len(rsdtypes_str.split(' '))
    if len(rsd.shape)==2:
        nrows, ncols = rsd.shape
        assert ncols == ncols_from_string, 'ncols from rsdtypes_str and rsd ncols do not match'
    else:
        nrows = rsd.shape[0]
        ncols = ncols_from_string

    lines = rsdheader.split('\n')
    # add following lines to header string to conform to radialshort data type
    lines.append('%' + 'TableType: LLUV RDL7')
    lines.append('%' + 'TableColumns: %d' % ncols)
    lines.append('%' + 'TableColumnTypes: %s' % rsdtypes_str)
    lines.append('%' + 'TableRows: %d' % nrows)
    lines.append('%TableStart:')
    lines.append('%%   Longitude   Latitude    U comp   V comp  VectorFlag    Spatial     Velocity    '+\
                 'Velocity  Velocity Spatial  X Distance  Y Distance   Range   Bearing   Velocity  '+\
                 'Direction   Spectra')
    lines.append('%%     (deg)       (deg)     (cm/s)   (cm/s)  (GridCode)    Quality     Maximum     '+\
                 'Minimum    Count    Count      (km)        (km)       (km)    (True)    (cm/s)     '+\
                 '(True)    RngCell')
    return '\n'.join(lines)

def unique_rows(a):
    # http://stackoverflow.com/questions/8560440/removing-duplicate-columns-and-rows-from-a-numpy-2d-array?lq=1
    a = numpy.ascontiguousarray(a)
    unique_a = numpy.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def cell_intersect(rngbear1, rngbear2):
    """ Return rows that match range and bearing data between the two nx2 matrices.

    Parameters
    ----------
    rngbear1 : nx2 array
       The first array of range and bearings that define each cell. 
       The first column is range cells.  The second column are the bearings.
    rngbear2 : nx2 array
       The second array of range and bearings that define each cell. 
       The first column is range cells.  The second column are the bearings.
    
    Return
    ------
    (rows1, rows2) : tuple of nx1 arrays 
       The rows where range and bearing cell are the same between two input matrices.
    """
    rows1 = []; rows2=[]
    for irow, cell in enumerate(rngbear1):
        rngcell, bearing = cell
        xrow = numpy.where( (rngbear2[:, 0] == rngcell) & \
                            (rngbear2[:, 1] == bearing) )[0]
        if xrow.size == 0:
            continue
        rows1.append(irow); rows2.append(xrow)
    rows1 = numpy.squeeze(rows1)
    rows2 = numpy.squeeze(rows2)
    return (rows1, rows2)


def compass2uv(wmag, wdir):
    """ Vector conversion from mag and direction (wmag,wdir) to x,y
    vector components (u,v)

    If inputs are lists, it is cast into an arrays.

    Parameters
    ----------
    wmag : array-like, same size as wdir
       The magnitude of the vector. 
    wdir : array-like, same size as wmag
       The compass direction, Clockwise from y-axis or North equals 0/360 deg

    Returns
    -------
    (u,v) : tuple of array-like u and v vectors the same size and shape as inputs.
       The x,y vector components. 

    >>> compass2uv(1.0, 0.0)
    (0.0, 1.0)
    >>> compass2uv([1., 1., 1., 1.], [0., 90., 180., 270.])
    (array([ 0.,  1.,  0., -1.]), array([ 1.,  0., -1., -0.]))
    
    """
    # calculate horizontal vector components (u,v) from magnitude and compass direction
    # cast the inputs into numpy.array
    wmag = numpy.array(wmag)
    wdir = numpy.array(wdir)
    assert wmag.shape == wdir.shape, 'wmag and wdir must be same size and shape'
    
    r = numpy.pi/180.
    u = wmag*numpy.sin(wdir*r)
    v = wmag*numpy.cos(wdir*r)
    return (u,v)



# for testing
if __name__ == '__main__':
    # 
    # datadir = sys.argv[1]
    # patterntype = sys.argv[2]
    # patterntype = 'MeasPattern' 
    # patterntype = 'IdealPattern'
    ifn = os.path.join('.', 'test', 'files', 'codar_raw', \
                   'Radialmetric_HATY_2013_11_05', \
                   'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    xd, xtypes_str = weighted_velocities(d, types_str)
    rsd, rsdtypes_str = generate_radialshort_array(d, types_str)
