#!/usr/bin/env python
# 
# Last modified: Time-stamp: <2017-07-06 20:17:04 codar>
""" CODAR Utilities 

"""
import sys
import os
import re
import fnmatch
import datetime

import geopy
import geopy.distance

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

def generate_radialshort_array(xd, xtypes_str, header, table_type='LLUV RDL7'):
    """Generates radialshort (rsd) data array.

    This function generates radialshort data (rsd) array based on data
    from weighted average data (xd).  It is the output CSV data table
    (middle) of CODAR LLUV format with header, middle, and footer. 

    If xd is empty, return empty rsd for writing out empty LLUV files.

    Each LATD, LOND is computed using Vincenty's algorithm for
    destination along a BEAR (NCW) and a distance of RNGE (km) from
    point of origin (lat1,lon1).  Vincenty's GC is great circle
    distance on an WGS-84 ellipsoid model between two points. This
    requires having the CODAR Site location, and range resolution from
    the header data or config data.
    
    Parameters
    ----------
    xd : ndarray
       QC'd and weighted average radials for each range, bearing where data were found
    xtypes_str : string 
        The order and key-labels for each column of xd array.
    header : string
        The header for the current date/time radialmetric required for site origin and range resolution
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

    if xd.size == 0:
        return numpy.array([]), rsdtypes_str

    # read header that match '%(k): (v)\n' pairs on each line
    m = re.findall(r'^(%.*):\s*(.*)$', header, re.MULTILINE)
    for k,v in m:
        ### print k+', '+v
        if k == '%TimeStamp':
            #sample_dt = scanf_datetime(v, fmt='%Y %m %d %H %M %S')
            pass
        elif k == '%Origin':
            lat1, lon1 = [float(x) for x in v.split()]
        elif k == '%RangeResolutionKMeters':
            range_resolution = float(v)
        elif k == '%TableStart':
            break

    # order of columns and labels for output data
    rsc = get_columns(rsdtypes_str)
    nrows,_ = xd.shape
    ncols = len(rsc)
    # Initialize new array for radial shorts
    rsd = numpy.ones(shape=(nrows,ncols))*numpy.nan
    rscol = numpy.array([rsc['VFLG'], rsc['SPRC'], rsc['BEAR'], rsc['VELO'], rsc['ESPC'], rsc['MAXV'], rsc['MINV'], rsc['EDVC'], rsc['ERSC']])

    xc = get_columns(xtypes_str)
    xcol = numpy.array([xc['VFLG'], xc['SPRC'], xc['BEAR'], xc['VELO'], xc['ESPC'], xc['MAXV'], xc['MINV'], xc['EDVC'], xc['ERSC']])

    # deal xd data into rsd by columns
    rsd[:,rscol] = xd[:,xcol]

    ############################
    # computations for filling in other columns of radialshort data 
    ############################

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
    rnge = range_resolution * rsd[:,rsc['SPRC']]
    rsd[:,rsc['RNGE']]=rnge
    #
    # Vincenty Great Circle destination point (LATD, LOND) based on rnge, bear from site origin
    origin = geopy.Point(lat1,lon1)
    pts = numpy.array([geopy.distance.vincenty(kilometers=r).destination(origin, b)[0:2] for (r,b) in zip(rnge,bear)])
    latd, lond = pts[:,0], pts[:,1]

    rsd[:,rsc['LATD']]=latd
    rsd[:,rsc['LOND']]=lond

    (xdist, ydist) = compass2uv(rnge,bear)
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

def run_LLUVMerger(ifn, ofn):
    """ Run CODAR's LLUVMerger app in subprocess """

    import subprocess

    # options mimic hourly radial merge for 
    cmdstr = '/Codar/SeaSonde/Apps/Bin/LLUVMerger -span=2.5 -lluvtype=i -angres=5 -angalign=2 -angmethod=short -method=average -minvect=2 -velcount -diag=4 -source='+ifn+' -output='+ofn

    print cmdstr
    subprocess.call(cmdstr, shell=True)
    return


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
