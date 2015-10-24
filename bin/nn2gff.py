#! /bin/env python3

'''
********************************************************************************
    (c) Andrew Robinson (andrew.robinson@latrobe.edu.au) 2013
        La Trobe University &
        Life Sciences Computation Centre (LSCC, part of VLSCI)
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Library General Public License (LGPL)
    as published by the Free Software Foundation; either version 2 of
    the License, or (at your option) any later version.

    ScienceScripts is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Library General Public License for more details.

    You should have recieved a copy of the GNU Library General Public
    License along with ScienceScripts; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
    USA
********************************************************************************

Created pm 08/02/2013

@author: arobinson
'''

from __future__ import print_function
from __future__ import division
from functools import reduce
import sys
import getopt
from Bio import SeqIO

def main(argv):

    # argument defaults
    minSize = 10

    delimiter = "\t"
    outrec = [
        '',        # 0, seq id
        'nn2gff',  # 1, program name
        'region',  # 2, feature type
        '0',       # 3, start idx (1-based)
        '0',       # 4, end idx (1-based inc)
        '.',       # 5, score
        '.',       # 6, strand
        '.',       # 7. phase
        '',        # 8. attributes
        ]

    # parse arguments
    try:
        opts, args = getopt.getopt(argv, 'hm:', [])
    except getopt.GetoptError:
        sys.stderr.write('nn2gff.py -m <min n count> <inputfile.fa>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print((
                "A simple script that searches FastA sequences for strings of "
                "N's. It outputs strings in GFF3 format on the standard "
                "output stream.\n\n"
                "nn2gff.py -m <min n count> <inputfile.fa>\n\n"
                "-m -- the minimum number of N's to consider a string\n"
                "<inputfile.fa> -- a text filename containing FastA formatted "
                "sequences\n"
                ))
            sys.exit()
        elif opt in ("-m",):
            minSize = int(arg)
    if len(args) == 1:
        filename = args[0]
    else:
        sys.stderr.write('Error: please specify only one file\n')
        sys.exit(3)

    # Do the work
    with open(filename, 'rU') as handle:
        # header
        print("##gff-version 3")

        # for each sequence
        for record in SeqIO.parse(handle, 'fasta'):
            stream = False  # True when we are in a string of N's
            streamstart = 0
            baseidx = 1
            outrec[0] = record.id
            for base in record.seq:
                if stream:  # in a string of N's
                    if base.upper() != "N":
                        stream = False
                        blocklen = baseidx - streamstart
                        if blocklen >= minSize:
                            outrec[3] = str(streamstart)
                            outrec[4] = str(baseidx)
                            outrec[8] = "length={}".format(blocklen)
                            print(delimiter.join(outrec))
                else:  # outside a string of N's
                    if base.upper() == 'N':
                        stream = True
                        streamstart = baseidx

                baseidx += 1

if __name__ == "__main__":
    main(sys.argv[1:])
