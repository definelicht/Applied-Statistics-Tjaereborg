from __future__ import division
import math, numpy as np, read_csv, root_utilities
from ROOT import TGraph

################################################################################

class Measurement:

  def __init__(self, parent, data):
    self.parent = parent
    self.raw = np.array(data,dtype="float")
    # Gaussian
    self.mean = np.mean(data)
    self.std = np.std(data)
    # Fourier
    self.n = data.shape[0]
    self.fourier_y = np.fft.fft(data,parent.next_power>>1) / self.n

def graph_rawdata(obj, line_color=4):
  graph = TGraph(obj.n,obj.parent.time,obj.raw)
  graph.SetLineColor(line_color)
  return graph

################################################################################

class WindData:

  def __init__(self, filename):
    # Read header info
    header = read_csv.Header(filename)
    self.height = []
    for h in header:
      try:
        self.height.append(float(h[2:len(h)-1]))
      except:
        continue
    self.height = np.array(self.height,dtype="float")
    self.n_heights = len(self.height)
    # Read data
    data = read_csv.ReadCsv(filename,skip=True)
    self.time = np.array(data[:,0],dtype="float")
    self.n = len(self.time)
    data = data[:,1:]
    self.t_lim = [np.min(self.time),np.max(self.time)]
    self.wind_lim = [np.min(data),np.max(data)]
    self.next_power = 2**(int(math.ceil(math.log(self.n, 2))))
    self.fourier_x = 0.5*35.0*np.linspace(0,1,self.next_power>>1)
    self.measurements = []
    for i in range(self.n_heights):
      self.measurements.append(Measurement(self,data[:,i]))

################################################################################

folder = "/Users/johannes/Dropbox/Applied Statistics/Project 2/Data/"

anemometer = WindData(folder + "windspeed.data")
sonicmeter = WindData(folder + "sonic.data")

################################################################################