#!/usr/bin/env python
import argparse
import json
import numpy
import pandas
import seaborn
from matplotlib import pyplot

seaborn.set()
seaborn.set_style("white")
seaborn.despine()

scpb = seaborn.color_palette("Blues")
scpr = seaborn.color_palette("Reds")
scpg = seaborn.cubehelix_palette(rot=-.4)

base_resolution = 100000

class Plotter:
  def __init__(self, mode, path, resolution):
    self.mode = mode
    self.path = path
    self.resolution = resolution
  # __init__

  def ParseJson(self):
    return json.load(open(self.path[0]))
  # ParseJson

  def DrawPile(self):
    try:
      data = self.ParseJson()
    except Exception:
      print("[tarantula::] error: wrong format")
      exit()

    for id in data:
      pile = data[id]

      figure, ax = pyplot.subplots(1, 1, figsize = (25, 5))

      ax.plot(range(len(pile["data"])), pile["data"], label = "data", color = scpb[2])
      ax.set_title(pile["id"])
      figure.text(0.5, 0.05, "base", ha = "center")
      figure.text(0.05, 0.5, "coverage", va = "center", rotation = "vertical")
      pyplot.legend(loc = "best")
      pyplot.savefig(str(pile["id"]) + ".png", format = 'png')
      pyplot.close(figure)
  # DrawPile

  def ParsePaf(self):
    paf1 = open(self.path[0], 'r')
    paf2 = open(self.path[1], 'r')
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

  def ParseSam(self):
    sam = open(self.path[0], 'r')
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

  def ParseNumpy(self):
    return numpy.load(self.path[0])
  # ParseNumpy

  def DrawHeatmap(self):
    if (self.resolution < base_resolution):
      return

    try:
      matrix = self.ParseNumpy()
    except Exception:
      try:
        matrix = self.ParseSam()
      except Exception:
        try:
          matrix = self.ParsePaf()
        except:
          print("[tarantula::] error: wrong format")
          exit()

    shrink_size = self.resolution // base_resolution
    heatmap_size = matrix.shape[0] // shrink_size
    deleted = [matrix.shape[0] - 1 - i for i in range(matrix.shape[0] % shrink_size)]

    heatmap = numpy.delete(matrix, deleted, axis = 0)
    heatmap = numpy.delete(heatmap, deleted, axis = 1)
    heatmap = heatmap.reshape((heatmap_size, shrink_size, heatmap_size, shrink_size)).sum(axis = 1).sum(axis = 2)
    heatmap = numpy.clip(heatmap, 0, 1000)

    fig, ax = pyplot.subplots(figsize = (80, 60))

    ax = seaborn.heatmap(heatmap,
      xticklabels = False,
      yticklabels = False,
      cmap = seaborn.cm.rocket_r,
      ax = ax)

    pyplot.savefig("heatmap_" + str(self.resolution) + ".png", format = 'png')
    pyplot.close()
  # DrawHeatmap

  def Run(self):
    if (self.mode == "heatmap"):
      self.DrawHeatmap()
    elif (self.mode == "pile"):
      self.DrawPile()
    return
  # Run
# Plotter

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description = "Plotter is a tool for drawing heatmaps and pile-o-grams",
    formatter_class = argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("mode",
    help = "draw either the [heatmap] or [pile]-o-grams")
  parser.add_argument("path",
    help = "heatmap: SAM | 2x PAF | npy - pile: JSON",
    nargs = "*")
  parser.add_argument("--resolution",
    help = "heatmap resolution in bp",
    default = base_resolution)

  args = parser.parse_args()
  plotter = Plotter(args.mode, args.path, args.resolution)
  plotter.Run()
# __main__
