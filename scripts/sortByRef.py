#!/usr/bin/env python

#Copyright (c) 2013, Regents of the University of California
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#1. Redistributions of source code must retain the above copyright notice,
#this list of conditions and the following disclaimer.
#
#2. Redistributions in binary form must reproduce the above copyright notice,
#this list of conditions and the following disclaimer in the documentation
#and/or other materials provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os
import sys
import argparse

def main():
    """ 
    Sorts lines of the input file INFILE according
    to the reference contig order specified by the
    reference dictionary REF_DICT (.fai file).
    The sort is stable. If --pos option is not specified,
    it is assumed that the contig name is the first
    field in each line.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pos', action='store', type=int, default=0)
    parser.add_argument('--tmp_dir', action='store', type=str, default='/tmp')
    parser.add_argument('input_file', action='store')
    parser.add_argument('output_file', action='store')
    parser.add_argument('ref_dict', action='store')
    args = parser.parse_args()

    # using dict and index is hacky 
    ref_dict = {} 
    index = 0 
    with open(args.ref_dict) as ref_file:
        for line in ref_file:
            contig = line.split('\t')[0] 
            ref_dict[contig] = index
            index = index + 1
    if args.input_file == '-':
        temp_outputs = {}
        in_file = sys.stdin 
        for line in in_file:
            fields = line.split('\t')
            contig = fields[args.pos].split(':')[0]
            if contig in ref_dict:
                order = ref_dict[contig]
            else:
                ref_dict[contig] = index
                order = index
                index = index + 1
            if order in temp_outputs: 
                tmp_file_handle = temp_outputs[order]
            else:
                tmp_name = args.tmp_dir + '/sortByRef_{}.tmp'.format(order)
                tmp_file_handle = open(tmp_name, 'w') 
                temp_outputs[order] = tmp_file_handle
            tmp_file_handle.write(line)
    else:
        temp_outputs = {}
        with open(args.input_file) as in_file:
            for line in in_file:
                fields = line.split('\t')
                contig = fields[args.pos].split(':')[0]
                if contig in ref_dict:
                    order = ref_dict[contig]
                else:
                    ref_dict[contig] = index
                    order = index
                    index = index + 1
                if order in temp_outputs: 
                    tmp_file_handle = temp_outputs[order]
                else:
                    tmp_name = args.tmp_dir + '/sortByRef_{}.tmp'.format(order)
                    tmp_file_handle = open(tmp_name, 'w') 
                    temp_outputs[order] = tmp_file_handle
                tmp_file_handle.write(line)

    for file_handle in temp_outputs.values():
        file_handle.close()

    if args.output_file == '-':
        out_file = sys.stdout 
        for order_num in range(index):
            if order_num in temp_outputs: 
                tmp_name = args.tmp_dir + '/sortByRef_{}.tmp'.format(order_num)
                with open(tmp_name) as tmp_file:
                    out_file.write(tmp_file.read())
                os.remove(tmp_name)                    
            else:
                continue
    else: 
        with open(args.output_file, 'w') as out_file:
            for order_num in range(index):
                if order_num in temp_outputs: 
                    tmp_name = args.tmp_dir + '/sortByRef_{}.tmp'.format(order_num)
                    with open(tmp_name) as tmp_file:
                        out_file.write(tmp_file.read())
                    os.remove(tmp_name)                    
                else:
                    continue


if __name__ == '__main__':
    main()
