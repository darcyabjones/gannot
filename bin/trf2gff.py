#! /usr/bin/env python3

from __future__ import print_function

program = "trf2gff"
version = "0.1.0"
author = "Darcy Jones"
date = "1 October 2015"
email = "darcy.ab.jones@gmail.com"
short_blurb = (
    'A simple script that parses output from '
    '[tandem repeat finder (TRF)](https://tandem.bu.edu/trf/trf.html) '
    'and returns the corresponding GFF3 file.\n'
    )
license = (
    '{program}-{version}\n'
    '{short_blurb}\n\n'
    'Copyright (C) {date},  {author}'
    '\n\n'
    'This program is free software: you can redistribute it and/or modify '
    'it under the terms of the GNU General Public License as published by '
    'the Free Software Foundation, either version 3 of the License, or '
    '(at your option) any later version.'
    '\n\n'
    'This program is distributed in the hope that it will be useful, '
    'but WITHOUT ANY WARRANTY; without even the implied warranty of '
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the '
    'GNU General Public License for more details.'
    '\n\n'
    'You should have received a copy of the GNU General Public License '
    'along with this program. If not, see <http://www.gnu.org/licenses/>.'
    )

license = license.format(**locals())

############################ Import all modules ##############################

import os
import re
import argparse
import sys
from collections import defaultdict

################################## Classes ###################################

class Feature(object):
    """ . """
    sep = '\t'
    attr_sep = ';'
    undefined = '.'
    column_order = [
        'seqid',
        'source',
        'type',
        'start',
        'end',
        'score',
        'strand',
        'phase',
        'attributes'
        ]

    def __init__(self, **kwargs):
        self.values = dict(zip(self.column_order, [None] * len(self.column_order)))
        self.values.update(kwargs)
        return

    def __str__(self):
        column_strings = []
        for col in self.column_order[:-1]:
            val = self.values[col]
            if col == 'strand' and val is not None:
                val = '-' if val == -1 else '+'
            column_strings.append(str(val) if val is not None else self.undefined)
        attrs = {k: v for k, v in self.values.items() if k not in self.column_order}
        column_strings.append(self._attributes(attrs))
        return self.sep.join(column_strings)

    def __repr__(self):
        return 'Feature({})'.format(', '.join(
                [self._attribute(k, v) for k, v in self.values.items()]
                ))

    def __getitem__(self, key):
        return self.values[key]

    def __setitem__(self, key, val):
        self.values[key] = val

    def _attribute(self, key, val):
        if type(val) not in {list, set, tuple}:
            val = [val]
        val = [str(v) for v in val]
        return '{}={}'.format(key, ','.join(val))

    def _attributes(self, attrs):
        attrs = [self._attribute(k, v) for k, v in attrs.items()]
        return self.attr_sep.join(attrs)


################################# Functions ##################################
fasta_name_line = re.compile(r">")
seqid_regex = re.compile(r"@|Sequence:", flags=re.IGNORECASE)
split_regex = re.compile(r"\s+")
data_regex = re.compile(
    r"+\s+".join([
        r'\d',
        r'\d',
        r'\d',
        r'[\d.]',
        r'\d',
        r'\d',
        r'\d',
        r'\d',
        r'\d',
        r'\d',
        r'\d',
        r'\d',
        r'[\d.]',
        r'\w',
        r'\w',
        ]) + r'+'
    )

def parse_fasta(handle):
    """ . """
    records = list()
    record = dict()
    for line in handle:
        if fasta_name_line.match(line) is not None:
            if len(record) > 0:
                records.append(record)
                record = dict()
            name = line.strip().lstrip('>').strip()
            record['id'] = name.split()[0]
            record['sequence'] = ''
        elif line.split() != "":
            record['sequence'] += line.strip().upper()
    records.append(record)
    return records

def reverse_complement(string):
    map_ = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        }
    complement = str()
    for base in reversed(string):
        complement += map_[base.upper()]
    return complement

def parse_trf(handle):
    """ . """
    columns = [
        ('start', int),
        ('end', int),
        ('period_size', int),
        ('copy_number', float),
        ('consensus_size', int),
        ('pmatch', int),
        ('pindel', int),
        ('score', int),
        ('A', int),
        ('C', int),
        ('G', int),
        ('T', int),
        ('entropy', float),
        ('consensus', str),
        ('repeat_sequence', str),
        ]
    seqid = None
    for line in handle:
        if seqid_regex.match(line):
            seqid = seqid_regex.sub('', line).strip()
            seqid = split_regex.split(seqid)[0]
        elif data_regex.match(line):
            line = split_regex.split(line.strip())
            record = {'seqid': seqid}
            for val, (col, type_) in zip(line, columns):
                record[col] = type_(val)
            yield record



################################# main code ##################################

def main(infile, outfile, repeat):
    def in_():
        if infile == '-':
            return sys.stdin
        else:
            return open(infile, 'rU')

    def out_():
        if outfile == '-':
            return sys.stdout
        else:
            return open(outfile, 'w')

    if repeat is not None:
        with open(repeat, 'r') as handle:
            sequences = parse_fasta(handle)
        repeats = {s['sequence']: s['id'] for s in sequences}
        repeats_rc = {reverse_complement(s['sequence']): s['id'] for s in sequences}

    seen_repeats = defaultdict(int)
    with in_() as inhandle, out_() as outhandle:
        for line in parse_trf(inhandle):
            strand = None
            if repeat is not None:
                fwd = line['consensus'] in repeats
                rev = line['consensus'] in repeats_rc
                if not (fwd or rev):
                    continue
                elif fwd:
                    line['type'] = repeats[line['consensus']]
                    strand = 1
                elif rev:
                    line['type'] = repeats_rc[line['consensus']]
                    strand = -1
            if 'type' not in line:
                line['type'] = 'SSR'
            seen_repeats[line['consensus']] += 1
            line['ID'] = line['consensus'] + '-' + str(seen_repeats[line['consensus']])
            line['source'] = 'TRF'
            line['strand'] = strand
            feature = Feature(
                seqid=line['seqid'],
                start=line['start'],
                end=line['end'],
                score=line['score'],
                phase=None,
                Name=line['consensus'],
                **{k: v for k, v in line.items() if \
                    k not in Feature.column_order}
                )
            outhandle.write(str(Feature(**line)) + "\n")

############################ Argument Handling ###############################

if __name__== '__main__':
    arg_parser = argparse.ArgumentParser(
      description=license,
      )
    arg_parser.add_argument(
        "-i", "--infile",
        default='-',
        help="Input TRF dat file. Default is '-' (stdin)"
        )
    arg_parser.add_argument(
        "-o","--outfile",
        default='-',
        help="Output path to write GFF to. Default is '-' (stdout)"
        )
    arg_parser.add_argument(
        "-r", "--repeat",
        default=None,
        help=("Fasta formatted sequences of repeats to include in the GFF. "
              "Fasta ids will be used as feature type. Default is to "
              "write all repeats with the feature type as 'SSR'.")
        )
    args = arg_parser.parse_args()

    main(**args.__dict__)
