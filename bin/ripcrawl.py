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

def main(argv, verbose=True):
    # argument defaults
    windowSize = 500
    incrementSize = 25
    headerRow = False
    maxNCount = -1

    delimiter = "\t"
    naValue = 'N/a'
    decimalPlaces = 5

    formatString = "{0:." + str(decimalPlaces) + "f}"
    header_row = [
        'seqid',
        'start',
        'end',
        'seq_length',
        'n_count',
        'product_index',
        'substrate_index',
        'composite_RIP_index',
        'gc_content'
        ]
    rows = list()
    if verbose:
        def printer(cols):
            sys.stdout.write(delimiter.join(cols) + '\n')
            rows.append(dict(zip(header_row, cols)))
            return
    else:
        def printer(cols):
            rows.append(dict(zip(header_row, cols)))
            return

    # parse arguments
    try:
        opts, args = getopt.getopt(argv, "htw:i:n:", [])
    except getopt.GetoptError:
        print('ripcrawl.py [-t] -n <max n count> -w <window size> -i <increment size> <inputfile.fa>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('ripcrawl.py [-t] -n <max n count> -w <window size> -i <increment size> <inputfile.fa>')
            sys.exit()
        elif opt in ("-i",):
            incrementSize = int(arg)
        elif opt in ("-w",):
            windowSize = int(arg)
        elif opt == "-t":
            headerRow = True
        elif opt == "-n":
            maxNCount = int(arg)

    if maxNCount < 1:
        maxNCount = int(windowSize * 0.1)

    if len(args) == 1:
        filename = args[0]
    else:
        sys.stderr.write("Error: please specify only one file\n")
        sys.exit(3)

    # Check windowSize and increment
    if windowSize / incrementSize != windowSize // incrementSize:
        sys.stderr.write("Error: window size must be a multiple of increment size\n")
        sys.exit(4)

    bucketCount = windowSize / incrementSize

    # Do the work
    if headerRow:
        printer(header_row)
    # for each sequence
    for record in SeqIO.parse(filename, 'fasta'):
        pairTotals = {
            'AA': 0.0, 'AT': 0.0, 'AC': 0.0, 'AG': 0.0,
            'TA': 0.0, 'TT': 0.0, 'TC': 0.0, 'TG': 0.0,
            'CA': 0.0, 'CT': 0.0, 'CC': 0.0, 'CG': 0.0,
            'GA': 0.0, 'GT': 0.0, 'GC': 0.0, 'GG': 0.0,
            }
        baseTotals = {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0,}
        incrementPairTotals = list()
        incrementBaseTotals = list()
        lastbase = 'N'
        baseidx = 0
        seqlen = len(record.seq)

        # for each base in sequence
        for base in record.seq:
            base = base.upper()
            # sum the new pair
            if base != 'N':
                baseTotals[base] += 1
                if lastbase != 'N':
                    pairTotals[lastbase + base] += 1
            # end of increment?
            if baseidx % incrementSize == 0 and baseidx != 0:
                incrementPairTotals.append(pairTotals)
                incrementBaseTotals.append(baseTotals)

                # Ready to output a full window?
                if len(incrementPairTotals) == bucketCount:
                    # Sum buckets
                    def l(x, y): return dict((k, v + y[k]) for k, v in x.items())
                    windowPairTotals = reduce(l, incrementPairTotals)
                    windowBaseTotals = reduce(l, incrementBaseTotals)

                    # calculate stats
                    productIndex = naValue
                    substrateIndex = naValue
                    compositeRIPIndex = naValue
                    gcContent = naValue

                    # ncount
                    nCount = windowSize - sum(windowPairTotals.values())

                    # product index
                    if nCount <= maxNCount:
                        if windowPairTotals['AT'] > 0:
                            productIndexF = windowPairTotals['TA'] / windowPairTotals['AT']
                            productIndex = formatString.format(productIndexF)

                        # substrate index
                        AC_GT = windowPairTotals['AC'] + windowPairTotals['GT']
                        if AC_GT > 0:
                            substrateIndex = (windowPairTotals['CA'] + windowPairTotals['TG']) / AC_GT

                            # composite rip index
                            if windowPairTotals['AT'] > 0:
                                compositeRIPIndex = formatString.format(productIndexF - substrateIndex)

                            substrateIndex = formatString.format(substrateIndex)

                        # GC content
                        ACGT = windowBaseTotals['A'] + windowBaseTotals['T'] + windowBaseTotals['C'] + windowBaseTotals['G']
                        if ACGT > 0:
                            gcContent = formatString.format((windowBaseTotals['C'] + windowBaseTotals['G']) / ACGT)

                    # print the results
                    printer([
                        record.id,  # seq id
                        str(baseidx - windowSize),  # window Start
                        str(baseidx - 1),  # window end
                        str(seqlen),  # sequence length
                        str(nCount),  # total N's in window
                        str(productIndex),
                        str(substrateIndex),
                        str(compositeRIPIndex),
                        str(gcContent),
                        ])
                    # remove old bucket
                    incrementPairTotals = incrementPairTotals[1:]
                    incrementBaseTotals = incrementBaseTotals[1:]

                # reset counters
                pairTotals = {
                    'AA': 0.0, 'AT': 0.0, 'AC': 0.0, 'AG': 0.0,
                    'TA': 0.0, 'TT': 0.0, 'TC': 0.0, 'TG': 0.0,
                    'CA': 0.0, 'CT': 0.0, 'CC': 0.0, 'CG': 0.0,
                    'GA': 0.0, 'GT': 0.0, 'GC': 0.0, 'GG': 0.0,
                    }
                baseTotals = {'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0,}
            # update counters and trackers
            baseidx += 1
            lastbase = base
    return rows

if __name__ == "__main__":
    main(sys.argv[1:])
