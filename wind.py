from __future__ import division
import math, numpy as np, read_csv, root_utilities
from ROOT import TCanvas, TGraph, TGraphErrors, TF1, SetOwnership, TMultiGraph, TLegend

################################################################################

class Measurement:

  def __init__(self, parent, data, height):
    self.parent = parent
    self.height = height
    self.raw = np.array(data,dtype="float")
    # Gaussian
    self.mean = np.mean(data)
    self.std = np.std(data)
    # Fourier
    self.fourier_y = np.fft.fft(data,parent.next_power) / parent.n
    self.fourier_y = 2.0*(np.absolute(self.fourier_y[:parent.next_power>>1]))**2
    self.fourier_y = np.array(self.fourier_y,dtype="float")
    # Power
    power_fit = TF1("powerfit%i"%id(self),"[0]*x^[1]",
                    parent.power_cutoff[0],parent.power_cutoff[1])
    power_fit.SetParLimits(1,-5,5)
    power_fit.SetLineWidth(3)
    power_graph = TGraph(parent.n,parent.fourier_x,self.fourier_y)
    power_graph.Fit("powerfit%i"%id(self),"r")
    SetOwnership(power_graph,False)
    self.power_fit = power_fit

def graph_rawdata(obj, line_color=4):
  graph = TGraph(obj.n,obj.parent.time,obj.raw)
  graph.SetLineColor(line_color)
  return graph

def graph_power(obj,fitcolor=7,discardcolor=4):

  i = 0

  # Values before fit range
  before_x = []
  before_y = []
  while obj.parent.fourier_x[i] < obj.parent.power_cutoff[0]:
    before_x.append(obj.parent.fourier_x[i])
    before_y.append(obj.fourier_y[i])
    i += 1
  before_x = np.array(before_x,dtype="float")
  before_y = np.array(before_y,dtype="float")
  before_graph = TGraph(before_x.shape[0],before_x,before_y)
  SetOwnership(before_graph,False)
  before_graph.SetLineColor(discardcolor)

  # Values in fit range
  fit_x = []
  fit_y = []
  while obj.parent.fourier_x[i] <= obj.parent.power_cutoff[1]:
    fit_x.append(obj.parent.fourier_x[i])
    fit_y.append(obj.fourier_y[i])
    i += 1
  fit_x = np.array(fit_x,dtype="float")
  fit_y = np.array(fit_y,dtype="float")
  fit_graph = TGraph(fit_x.shape[0],fit_x,fit_y)
  SetOwnership(fit_graph,False)
  fit_graph.SetLineColor(fitcolor)

  # Values after fit range
  after_x = []
  after_y = []
  while i < obj.parent.n_fourier:
    after_x.append(obj.parent.fourier_x[i])
    after_y.append(obj.fourier_y[i])
    i += 1
  after_x = np.array(after_x,dtype="float")
  after_y = np.array(after_y,dtype="float")
  after_graph = TGraph(after_x.shape[0],after_x,after_y)
  SetOwnership(after_graph,False)
  after_graph.SetLineColor(discardcolor)

  # Add fit
  fit_graph.Fit("powerfit%i"%id(obj))

  # Create multigraph
  multigraph = TMultiGraph()
  multigraph.Add(before_graph)
  multigraph.Add(fit_graph)
  multigraph.Add(after_graph)

  return multigraph

################################################################################

class WindData:

  power_cutoff = [1.0e-1,3e0]

  def __init__(self, filename):
    # Read header info
    header = read_csv.Header(filename)
    self.heights = []
    for h in header:
      try:
        self.heights.append(float(h[2:len(h)-1]))
      except:
        continue
    self.heights = np.array(self.heights,dtype="float")
    self.n_heights = len(self.heights)
    # Read data
    data = read_csv.ReadCsv(filename,skip=True)
    self.time = np.array(data[:,0],dtype="float")
    self.n = len(self.time)
    data = data[:,1:]
    self.t_lim = [np.min(self.time),np.max(self.time)]
    self.wind_lim = [np.min(data),np.max(data)]
    self.next_power = 2**(int(math.ceil(math.log(self.n, 2))))
    self.fourier_x = 0.5*35.0*np.linspace(0,1,self.next_power>>1)
    self.n_fourier = len(self.fourier_x)
    self_fourier_lim = [np.min(self.fourier_x),np.max(self.fourier_x)]
    self.measurements = []
    for i in range(self.n_heights):
      self.measurements.append(Measurement(self,data[:,i],self.heights[i]))

################################################################################

folder = "/Users/johannes/Dropbox/Applied Statistics/Project 2/Data/"

anemometer = WindData(folder + "windspeed.data")
sonicmeter = WindData(folder + "sonic.data")

################################################################################

def draw_power_slopes(save=False):

  def get(measurements):
    n = len(measurements)
    heights = []
    slopes = []
    errors = []
    for m in measurements:
      heights.append(m.height)
      slopes.append(m.power_fit.GetParameter(1))
      errors.append(m.power_fit.GetParError(1))
    heights = np.array(heights,dtype="float")
    slopes = np.array(slopes,dtype="float")
    errors = np.array(errors,dtype="float")
    return (n,heights,slopes,errors)

  # Build the two graphs
  (ane_n,ane_heights,ane_slopes,ane_errors) = get(anemometer.measurements)
  (son_n,son_heights,son_slopes,son_errors) = get(sonicmeter.measurements)
  ane_graph = TGraphErrors(ane_n,ane_heights,ane_slopes,np.zeros(ane_n),
                           ane_errors)
  son_graph = TGraphErrors(son_n,son_heights,son_slopes,np.zeros(son_n),
                           son_errors)
  SetOwnership(ane_graph,False)
  SetOwnership(son_graph,False)

  # Aggregate analysis
  slopes = np.concatenate((ane_slopes,son_slopes))
  mean = [np.mean(ane_slopes),np.mean(son_slopes),np.mean(slopes)]
  std = [np.std(ane_slopes),np.std(son_slopes),np.std(slopes)]

  # Build multigraph
  canvas = TCanvas()
  SetOwnership(canvas,False)
  multigraph = TMultiGraph()
  SetOwnership(multigraph,False)
  multigraph.Add(ane_graph)
  multigraph.Add(son_graph)
  multigraph.SetTitle(";Height [m];Slope of power fit [1]")
  multigraph.Draw("AP")
  multigraph.GetYaxis().SetRangeUser(mean[2]-1,mean[2]+1)
  (x1,x2) = (ane_heights[0]-abs(ane_heights[0]),ane_heights[-1]+abs(ane_heights[-1]))
  lines = [
    root_utilities.DrawLine(x1=x1,x2=x2,
      y=mean[0],color=2,width=2,style=2),
    root_utilities.DrawLine(x1=x1,x2=x2,
      y=mean[1],color=4,width=2,style=2),
    root_utilities.DrawLine(x1=x1,x2=x2,
      y=mean[2],color=6,width=2,style=2),
    root_utilities.DrawLine(x1=x1,x2=x2,
      y=-5/3,color=8,width=2,style=1)
  ]

  # Build legend
  line_legend = TLegend(0.5,0.1,0.9,0.9-0.5)
  line_legend.SetFillColor(0)
  line_legend.AddEntry(lines[0],"Anemometer mean slope","L")
  line_legend.AddEntry(lines[1],"Sonicmeter mean slope","L")
  line_legend.AddEntry(lines[2],"Total mean slope","L")
  line_legend.AddEntry(lines[3],"Theoretical slope","L")
  line_legend.SetLineColor(1)
  line_legend.Draw()
  point_legend = TLegend(0.5,0.7,0.9,0.9)
  point_legend.SetFillColor(0)
  point_legend.AddEntry(ane_graph,"Anemometer slopes","PE")
  point_legend.AddEntry(son_graph,"Sonicmeter slopes","PE")
  point_legend.Draw()

  # Cosmetics
  root_utilities.SetLine(ane_graph,color=2,width=2,style=1)
  root_utilities.SetLine(son_graph,color=4,width=2,style=1)
  root_utilities.SetMarker(ane_graph,color=2,size=2,style=4)
  root_utilities.SetMarker(son_graph,color=4,size=2,style=5)

  canvas.Update()

  if save:
    canvas.SaveAs("plots/slopes.pdf")
  else:
    raw_input("Displaying graph. Press enter to continue.")

def draw_power_fit(obj,save=False):
  canvas = TCanvas()
  canvas.SetLogx()
  canvas.SetLogy()
  graph = graph_power(obj)
  graph.Draw("AL")
  if save:
    canvas.SaveAs("plots/power.pdf")
  else:
    raw_input("Displaying graph. Press enter to continue.")
