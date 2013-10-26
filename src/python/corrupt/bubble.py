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
