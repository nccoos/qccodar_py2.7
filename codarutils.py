#!/usr/bin/env python
# 
# Last modified: Time-stamp: <2015-06-01 17:54:44 haines>
""" CODAR Utilities 

"""
import sys
import os
import re

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
        print 'No Radial Data in '+ ifn
        return numpy.array([]), '', header, footer

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
