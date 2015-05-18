#
# Last modified: Time-stamp: <>
"""
CODAR File Utilities 
"""
import sys
import os
import re
import glob

import numpy
from datetime import datetime
from time import strptime
from StringIO import StringIO
from collections import defaultdict

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

def write_empty_output(ofn, header, footer):
    "write header and footer only, since no radial data"
    f = open(ofn, 'w')
    f.write(header)
    f.write(footer)
    f.close()

def write_output(ofn, header, d, footer):
    "write header, radialmetric data, and footer "
    f = open(ofn, 'w')
    f.write(header)
    numpy.savetxt(f, d)
    f.write(footer)
    f.close()

def read_lluv_file(ifn):
    lines = load_data(ifn)

    m=re.match(r'(?P<header>(%.*\n)*)(?P<middle>([\d\s-].*\n)*)(?P<tail>(%.*\n)*)', ''.join(lines))
    header  = m.group('header')
    footer = m.group('tail')

    ####
    # print ifn
    # print 'header: %d, tail: %d' % (len(header), len(footer))
    # print footer

    if len(footer)<=0:
        print 'No Radial Data in '+ ifn
        return '', '', header, footer

    # read header that match '%(k): (v)\n' pairs on each line
    m = re.findall(r'^(%.*):\s*(.*)$', ''.join(lines), re.MULTILINE)
    for k,v in m:
        #### print k+', '+v
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
    
    ####
    # print types_str

    # identify column numbers for selected variables
    m2 = re.findall(r'\w{4}', types_str)
        
    # if len(footer)<=0:
    #     print 'No Radial Data in '+ ifn
    #     # empty array to append
    #     d = numpy.array([]).reshape(0,len(m2))
    #     return d, types_str, header, footer

    # read data from string of lines but make it behave like a file object with StringIO
    s = StringIO(''.join(lines))
    s.seek(0) # ensures start posn of file-like string s
    d = numpy.loadtxt(s, comments='%')
    # lat, lon, u, v = numpy.loadtxt(s, usecols=(0,1,2,3), comments='%', unpack=True)
    return d, types_str, header, footer

def get_columns(types_str):
    # use dict to store column label and it's column number
    c = defaultdict(int)
    column_labels = types_str.strip().split(' ')
    m = re.findall(r'\w{4}', types_str)
    for label in column_labels:
        c[label]=m.index(label) # c['VFLG']=4
    return c
