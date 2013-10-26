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


