#Copyright (c) 2013, Regents of the University of California
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from random import random, randint
import sys

def mutate_base(base):
  bases = ['A', 'C', 'G', 'T']
  reverse_bases = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
  index = (reverse_bases[base] + randint(1,3)) % 4
  return bases[index]

def mutate_read(read):
  local_length = sublength
  if random() < likelihood:
    char_list = list(read)
    mutate_index = 50
    while local_length != 0:
      char_list[mutate_index] = mutate_base(char_list[mutate_index])
      mutate_index += 1
      local_length -= 1
    return "".join(char_list)
  else:
    return read

if __name__ == '__main__':
  input_file_name = sys.argv[1]
  output_file_name = sys.argv[2]
  if len(sys.argv) > 3:
    likelihood = float(sys.argv[3])
  else:
    likelihood = 0.04
  if len(sys.argv) > 4:
    sublength = int(sys.argv[4])
  else:
    sublength = 1
  wf = open(output_file_name, 'w')
  rf = open(input_file_name)

  counter = 0
  for line in rf:
    if counter % 4 == 1:
      wf.write(mutate_read(line))
    else:
      wf.write(line)
    counter += 1
  print "Done!"


