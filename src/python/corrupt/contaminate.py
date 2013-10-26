from random import random, randint
import sys

if __name__ == '__main__':
  input_file_name = sys.argv[1]
  contaminate_file_name = sys.argv[2]
  output_file_name = sys.argv[3]
  if len(sys.argv) > 4:
    likelihood = float(sys.argv[4])
  else:
    likelihood = 0.04
  wf = open(output_file_name, 'w')
  rf = open(input_file_name)
  cf = open(contaminate_file_name)

  counter = 0
  contam = False
  for line in rf:
    if counter % 4 == 0:
      if random() < likelihood:
        contam = True
      else:
        contam = False
    if contam:
      if counter % 4 == 0:
        wf.write(line)
        cf.next()
      else:
        wf.write(cf.next())
    else:
      wf.write(line)
    counter += 1
  print "Done!"
