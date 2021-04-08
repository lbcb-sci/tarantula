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

def clip_matrix(infile,outfile, vmax, vmin, value):
  df = pd.read_csv(infile)
  df.drop(['Total','Unnamed: 0'], axis=1, inplace=True)
  lens = df.shape
  name = outfile[:-4]
  print(name+'_clip_' +str(value)+'.csv')

  print(lens)
  print(lens[0])
  for i in range(lens[0]):
    for r in range(lens[1]):
      if (df.iloc[i][r] > 1000):
        df.iat[i, r] = 1000
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
  name = outfile[:-4]
  print(name)
  df.to_csv(name+'_clip_' +str(value)+'.csv')

if __name__ == "__main__":
  arguments = len(sys.argv)-1
  longOptions = ['max','min','clip']
  options = "m:i:c:"
  max = -1
  min = -1
  clip = False
  try: 
    opts, args = getopt.getopt(sys.argv[1:], options, longOptions)
  except getopt.GetoptError:
    print('wrong params')
  for opt,arg in opts:
    if opt in ("-m", "--max"):
      max = int(arg)
    elif opt in ("-i", "--min"):
      min = int(arg)
    elif opt in ("-c", "--clip"):
      clip = True
      clip_value = int(arg)

  inFile = sys.argv[arguments-1]
  outFile = sys.argv[arguments]
  if (clip):
    clip_matrix(inFile, outFile, max, min, clip_value)
  else:
  
    plot(inFile, outFile, max, min)