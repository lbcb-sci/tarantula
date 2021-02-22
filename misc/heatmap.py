#!/usr/bin/env python
import numpy
import pandas
import seaborn
import sys
from matplotlib import pyplot

seaborn.set()
seaborn.set_style("white")
seaborn.despine()

base_resolution = 100000

def ParseSam(path):
  sam = open(path, 'r')
  references = {}
  line = ""
  while True: # parse header
    line = sam.readline()
    if (len(line) == 0 or line[0] != '@'):
      break
    if (line[1:3] == 'SQ'):
      reference_name = line.split('SN:')[1].split('\t')[0]
      reference_len = int(line.split('LN:')[1].split('\t')[0].rstrip())
      if (reference_len > 10000):
        references[reference_name] = reference_len

  references = dict(sorted(references.items(), key = lambda item: item[1], reverse = True))
  matrix_size = 0
  for reference_name, reference_len in references.items():
    references[reference_name] = matrix_size
    matrix_size += reference_len

  matrix_size = matrix_size // base_resolution + 1
  matrix = numpy.zeros(shape = (matrix_size, matrix_size))

  while True: # parse alignments
    _, flag, rhs_name, rhs_begin, quality, _, rhs_next_name, rhs_next_begin, *junk = line.split('\t')
    if (rhs_next_name == '='):
      rhs_next_name = rhs_name

    flag = int(flag)
    rhs_begin = int(rhs_begin)
    rhs_next_begin = int(rhs_next_begin)

    #   1 - multiple segments
    #   4 - segment unmapped
    #   8 - next segment unmapped
    #  64 - first segment
    # 256 - secondary alignment
    if (    (flag &   1) and \
        not (flag &   4) and \
        not (flag &   8) and \
            (flag &  64) and \
        not (flag & 256)):
      if (rhs_name      in references and \
          rhs_next_name in references):
        x = (rhs_begin      + references[rhs_name])      // base_resolution
        y = (rhs_next_begin + references[rhs_next_name]) // base_resolution
        matrix[x][y] += 1
        matrix[y][x] += 1

    line = sam.readline()
    if (len(line) == 0):
      break

  numpy.save("heatmap", matrix)

  return matrix
# ParseSam

def ParseNumpy(path):
  matrix = numpy.load(path)
  return matrix
# ParseNumpy

def Plot(matrix, resolution):
  if (resolution < base_resolution):
    return

  shrink_size = resolution // base_resolution
  heatmap_size = matrix.shape[0] // shrink_size
  deleted = [matrix.shape[0] - 1 - i for i in range(matrix.shape[0] % shrink_size)]

  heatmap = numpy.delete(matrix, deleted, axis = 0)
  heatmap = numpy.delete(heatmap, deleted, axis = 1)
  heatmap = heatmap.reshape((heatmap_size, shrink_size, heatmap_size, shrink_size)).sum(axis = 1).sum(axis = 2)
  heatmap = numpy.clip(heatmap, 0, 1000)

  ax = seaborn.heatmap(heatmap, \
    xticklabels = False, \
    yticklabels = False, \
    cmap = seaborn.cm.rocket_r)

  pyplot.savefig("heatmap_" + str(resolution) + ".png", format = 'png')
  pyplot.close()

if __name__ == "__main__":
  try:
    matrix = ParseNumpy(sys.argv[1])
  except Exception:
    try:
      matrix = ParseSam(sys.argv[1])
    except Exception:
      print("[tarantula::] error: wrong format")
      exit()

  for i in range(2, len(sys.argv)):
    Plot(matrix, int(sys.argv[i]))
