import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
import sys, getopt

def plot(infile, outfile, vmax, vmin):
  df = pd.read_csv(infile)
  df.drop(['Total','Unnamed: 0'], axis=1, inplace=True)
  fig, ax = plt.subplots(figsize = (80, 60))
  
  if (vmax != -1 and vmin != -1):
    ax = sb.heatmap(df, cmap = sb.cm.rocket_r, ax = ax, vmax = vmax, vmin = vmin)
  elif (vmax != -1):
    ax = sb.heatmap(df, cmap = sb.cm.rocket_r, ax = ax, vmax = vmax)
  elif (vmin != -1):
    ax = sb.heatmap(df, cmap = sb.cm.rocket_r, ax = ax, vmin = vmin)
  else:
    ax = sb.heatmap(df, cmap = sb.cm.rocket_r, ax = ax)
  plt.savefig(outfile, dpi=100)

if __name__ == "__main__":
  arguments = len(sys.argv)-1
  longOptions = ['max']
  options = "m:i:"
  max = -1
  min = -1
  try: 
    opts, args = getopt.getopt(sys.argv[1:], options, longOptions)
  except getopt.GetoptError:
    print('wrong params')
  for opt,arg in opts:
    if opt in ("-m", "--max"):
      max = int(arg)
    elif opt in ("-i", "--min"):
      min = int(arg)

  inFile = sys.argv[arguments-1]
  outFile = sys.argv[arguments]
  plot(inFile, outFile, max, min)