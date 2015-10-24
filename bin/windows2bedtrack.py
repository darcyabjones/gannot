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
from math import ceil
from math import floor

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
def strand_in(value):
    """ . """
    in_map = {
        '+': 1,
        '-': -1,
        '1': 1,
        '-1':-1
        }
    try:
        return in_map[value]
    except KeyError:
        raise ValueError('Invalid value "{}" in strand field.'.format(value))

def strand_out(value):
    """ . """
    out_map = {1: '+', -1: '-', '': ''}
    return out_map[value]

def field_list_in(value, type_=int):
    """ . """
    return [type_(v.strip()) for v in field.split(',')]

def field_list_out(value):
    """ . """
    return ','.join([str(v) for v in value])

def parse_bed(handle):
    """ . """
    columns = [
        ('seqid', str, str),
        ('start', int, str),
        ('end', int, str),
        ('name', str, str),
        ('score', float, str),
        ('strand', strand_in, strand_out),
        ('thick_start', int, str),
        ('thick_end', int, str),
        ('item_rgb', field_list_in, field_list_out),
        ('block_count', int, str),
        ('block_sizes', field_list_in, field_list_out),
        ('block_starts', field_list_in, field_list_out),
        ]
    na = ['N/a', '.', '-']

    for line in handle:
        line = line.strip().split('\t')
        record = dict()
        try:
            for val, (col, in_type, out_type) in zip(line, columns):
                record[col] = in_type(val)
            yield record
        except ValueError:
            if val in na:
                pass
            else:
                raise sys.exc_info()[1]

def write_bed(handle, record):
    """ . """
    columns = [
        ('seqid', str, str),
        ('start', int, str),
        ('end', int, str),
        ('name', str, str),
        ('score', float, str),
        ('strand', strand_in, strand_out),
        ('thick_start', int, str),
        ('thick_end', int, str),
        ('item_rgb', field_list_in, field_list_out),
        ('block_count', int, str),
        ('block_sizes', field_list_in, field_list_out),
        ('block_starts', field_list_in, field_list_out),
        ]
    sep = '\t'
    line = list()
    for col, in_type, out_type in columns:
        if col not in record:
            val = ''
        else:
            val = record[col]
        line.append(out_type(val))
    handle.write(sep.join(line).strip() + '\n')
    return

def write_bedgraph(handle, record):
    """ . """
    columns = [
        ('seqid', str, str),
        ('start', int, str),
        ('end', int, str),
        ('score', float, str),
        ]
    sep = '\t'
    line = list()
    for col, in_type, out_type in columns:
        if col not in record:
            val = ''
        else:
            val = record[col]
        line.append(out_type(val))
    handle.write(sep.join(line).strip() + '\n')
    return

def write_track_header(handle, record):
    """ . """
    line = ['track']
    for key, val in record.items():
        val = str(val)
        if ' ' in val:
            val = '"' + val + '"'
        line.append(str(key) + "=" + val)
    handle.write(' '.join(line) + '\n')


################################# main code ##################################

def main(infile, outfile, name):
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

    with in_() as in_handle, out_() as out_handle:
        if name:
            write_track_header(out_handle, {'name': name, 'type': 'bedGraph',})
        last_line = None
        for line in parse_bed(in_handle):
            if last_line is None:  # First record
                last_line = line
            else:
                interval = line['start'] - last_line['start']
                center = last_line['start'] + ceil(0.5 * (last_line['end'] - last_line['start']))
                last_line['start'] = center - floor(0.5 * interval)
                last_line['end'] = center + ceil(0.5 * interval) + 1
                write_bedgraph(out_handle, last_line)
                last_line = line
        center = last_line['start'] + 0.5 * (last_line['end'] - last_line['start'])
        last_line['start'] = center - 0.5 * interval
        last_line['end'] = center + 0.5 * interval
        write_bedgraph(out_handle, last_line)


############################ Argument Handling ###############################

if __name__== '__main__':
    arg_parser = argparse.ArgumentParser(
      description=license,
      )
    arg_parser.add_argument(
        "-i", "--infile",
        default='-',
        help="Input sorted bed window file. Default is '-' (stdin)"
        )
    arg_parser.add_argument(
        "-o","--outfile",
        default='-',
        help="Output path to write track to. Default is '-' (stdout)"
        )
    arg_parser.add_argument(
        "-n","--name",
        default=None,
        help="Name to call the track."
        )
    args = arg_parser.parse_args()

    main(**args.__dict__)
