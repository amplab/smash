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

import random
import sys

def bubble_bases(read):
  bases = ['A', 'C', 'G', 'T']
  bubble_insert = reduce(lambda x, y: x + bases[random.randint(0,3)], range(0, bubble_length), '')
  return read[0:50] + bubble_insert + read[50:bubble_index] + '\n'

def bubble_quality(line):
  return line[0:50] + reduce(lambda x, y: x + chr(random.randint(33, 126)), range(0, bubble_length), '') + line[50:bubble_index] + '\n'

if __name__ == '__main__':
  input_file_name = sys.argv[1]
  output_file_name = sys.argv[2]
  if len(sys.argv) > 3:
    likelihood = float(sys.argv[3])
  else:
    likelihood = 0.04
  if len(sys.argv) > 4:
    bubble_length = int(sys.argv[4])
  else:
    bubble_length = 1
  bubble_index = -(bubble_length + 1)

  wf = open(output_file_name, 'w')
  rf = open(input_file_name)

  header = rf.next()
  while header:
    try:
      corrupt = random.random() < likelihood
      wf.write(header)
      bases = rf.next()
      if corrupt:
        wf.write(bubble_bases(bases))
      else:
        wf.write(bases)
      plus = rf.next()
      wf.write(plus)
      quality = rf.next()
      if corrupt:
        wf.write(bubble_quality(quality))
      else:
        wf.write(quality)
      header = rf.next()
    except StopIteration:
      break
  rf.close()
  wf.close()
  print "Done!"
