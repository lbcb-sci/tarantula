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

def ParsePaf(path1, path2):
  paf1 = open(path1, 'r')
  paf2 = open(path2, 'r')
  references = {}
  connections = {}

  rp1 = {}
  line = ""
  while True:
    line = paf1.readline()
    if (len(line) == 0):
      break
    lhs_name, _, _, _, _, rhs_name, rhs_len, rhs_begin, *junk = line.split('\t')
    rhs_len = int(rhs_len)
    rhs_begin = int(rhs_begin)
    if (rhs_name not in references and rhs_len > 10000):
      references[rhs_name] = rhs_len
    if (lhs_name not in rp1):
      rp1[lhs_name] = (rhs_name, rhs_begin)

  rp2 = {}
  while True:
    line = paf2.readline()
    if (len(line) == 0):
      break
    lhs_name, _, _, _, _, rhs_name, rhs_len, rhs_begin, *junk = line.split('\t')
    rhs_len = int(rhs_len)
    rhs_begin = int(rhs_begin)
    if (rhs_name not in references and rhs_len > 10000):
      references[rhs_name] = rhs_len
    if (lhs_name not in rp2):
      rp2[lhs_name] = (rhs_name, rhs_begin)

  references = dict(sorted(references.items(), key = lambda item: item[1], reverse = True))
  matrix_size = 0
  for reference_name, reference_len in references.items():
    references[reference_name] = matrix_size
    matrix_size += reference_len

  matrix_size = matrix_size // base_resolution + 1
  matrix = numpy.zeros(shape = (matrix_size, matrix_size))

  for lhs_name, rhs in rp1.items():
    if (lhs_name not in rp2 or \
        rhs[0] not in references or \
        rp2[lhs_name][0] not in references):
      continue
    x = (rhs[1]           + references[rhs[0]])          // base_resolution
    y = (rp2[lhs_name][1] + references[rp2[lhs_name][0]]) // base_resolution
    matrix[x][y] += 1
    matrix[y][x] += 1

  numpy.save("heatmap", matrix)

  return matrix
# ParsePaf

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
  if (len(references) == 0):
    raise Error

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

  fig, ax = pyplot.subplots(figsize = (80, 60))

  ax = seaborn.heatmap(heatmap, \
    xticklabels = False, \
    yticklabels = False, \
    cmap = seaborn.cm.rocket_r, \
    ax = ax)

  pyplot.savefig("heatmap_" + str(resolution) + ".png", format = 'png')
  pyplot.close()
# Plot

if __name__ == "__main__":
  try:
    matrix = ParseNumpy(sys.argv[1])
  except Exception:
    try:
      matrix = ParseSam(sys.argv[1])
    except Exception:
      try:
        matrix = ParsePaf(sys.argv[1], sys.argv[2])
      except:
        print("[tarantula::] error: wrong format")
        exit()

  for i in range(2, len(sys.argv)):
    if (sys.argv[i].isdigit()):
      Plot(matrix, int(sys.argv[i]))
# __main__