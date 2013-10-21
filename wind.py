from __future__ import division
import math, numpy as np, read_csv, root_utilities
from ROOT import TCanvas, TGraph, TGraphErrors, TF1, SetOwnership, TMultiGraph, TLegend, TColor, TH1F

################################################################################

class Measurement:

  def __init__(self, parent, data, height):
    self.parent = parent
    self.height = height
    self.raw = np.array(data,dtype="float")
    # Gaussian
    self.mean = np.mean(data)
    self.std = np.std(data) / (parent.n-1)**0.5
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

  # Get data
  (ane_n,ane_heights,ane_slopes,ane_errors) = get(anemometer.measurements)
  (son_n,son_heights,son_slopes,son_errors) = get(sonicmeter.measurements)

  # Aggregate analysis
  mean = [np.mean(ane_slopes),np.mean(son_slopes),0]
  ane_error = 0
  son_error = 0
  for slope in ane_slopes: ane_error += (slope - mean[0])**2
  for slope in son_slopes: son_error += (slope - mean[1])**2
  ane_error = ane_error**0.5 / (ane_n-1)
  son_error = son_error**0.5 / (son_n-1)
  ane_errors = ane_error*np.ones(ane_n,dtype="float")
  son_errors = son_error*np.ones(son_n,dtype="float")
  (mean[2],error) = root_utilities.WeightedMean(
    np.concatenate((ane_slopes,son_slopes)),
    np.concatenate((ane_errors,son_errors))
  )

  # Build the two graphs
  ane_graph = TGraphErrors(ane_n,ane_heights,ane_slopes,np.zeros(ane_n),
                           ane_errors)
  son_graph = TGraphErrors(son_n,son_heights,son_slopes,np.zeros(son_n),
                           son_errors)
  SetOwnership(ane_graph,False)
  SetOwnership(son_graph,False)

  # Build multigraph
  canvas = TCanvas()
  SetOwnership(canvas,False)
  multigraph = TMultiGraph()
  SetOwnership(multigraph,False)
  multigraph.Add(ane_graph)
  multigraph.Add(son_graph)
  multigraph.SetTitle("Power fit slopes;Height [m];Slope of power fit [1]")
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
  line_legend.AddEntry(lines[0],"Anemometer: %.2f"%mean[0],"L")
  line_legend.AddEntry(lines[1],"Sonicmeter: %.2f"%mean[1],"L")
  line_legend.AddEntry(lines[2],"Total: %.2f"%mean[2],"L")
  line_legend.AddEntry(lines[3],"Theoretical: %.2f"%(-5/3),"L")
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

  raw_input("Displaying graph. Press enter to continue.")
  if save:
    canvas.SaveAs("plots/slopes.pdf")

################################################################################

def graph_rawdata(obj, line_color=4):
  graph = TGraph(obj.parent.n,obj.parent.time,obj.raw)
  SetOwnership(graph,False)
  graph.SetLineColor(line_color)
  return graph

################################################################################

def draw_rawdata(save=False):

  # Setup
  multigraph = TMultiGraph()
  multigraph.SetTitle("Raw wind speed data;Time [s];Wind [m/s]")
  SetOwnership(multigraph,False)
  color = 1
  legend = TLegend(0.8,0.7,0.95,0.95)

  # Add each measurement for the anemometer
  for measurement in anemometer.measurements:
    graph = graph_rawdata(measurement,line_color=color)
    multigraph.Add(graph)
    multigraph.Add(root_utilities.DrawLine(x1=anemometer.t_lim[0],
        x2=anemometer.t_lim[1],y=measurement.mean,color=color,style=2,width=2))
    legend.AddEntry(graph,"h = %.1fm"%measurement.height,"L")
    color += 1

  # Draw
  canvas = TCanvas()
  SetOwnership(canvas,False)
  multigraph.Draw("AL")
  legend.Draw()

  raw_input("Displaying graph. Press enter to continue.")
  if save: canvas.SaveAs("plots/rawdata.pdf")

################################################################################

def draw_power(obj, secondary=None, fitcolor=7, discardcolor=4, save=False):
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
  fit_graph.Fit("powerfit%i"%id(obj),"r")

  # Create multigraph
  multigraph = TMultiGraph()
  SetOwnership(multigraph,False)
  if secondary != None:
    secondary_graph = TGraph(obj.parent.n_fourier,obj.parent.fourier_x,
                             secondary.fourier_y)
    secondary_graph.SetLineColor(9)
    multigraph.Add(secondary_graph)
  multigraph.Add(before_graph)
  multigraph.Add(fit_graph)
  multigraph.Add(after_graph)

  # Draw
  canvas = TCanvas()
  SetOwnership(canvas,False)
  canvas.SetLogx()
  canvas.SetLogy()
  multigraph.SetTitle("Power spectrum for h = %.1fm;f [Hz];|Y(f)|^{2}"
                      % obj.height)
  multigraph.Draw("AL")

  # Legends
  legend = TLegend(0.6,0.6,0.88,0.88)
  legend.SetFillStyle(0)
  legend.SetLineColor(0)
  legend.AddEntry(before_graph,"Discarded data","L")
  legend.AddEntry(fit_graph,"Data used in fit","L")
  if secondary != None:
    legend.AddEntry(secondary_graph,"Data for h = %.1fm"%secondary.height,"L")
  legend.Draw()
  legend_prob = TLegend(0.12,0.12,0.7,0.3)
  root_utilities.StealthyLegend(legend_prob)
  legend_prob.AddEntry(obj.power_fit,"f(x) = %.2fx^{%.2f}" % 
                       (obj.power_fit.GetParameter(0),
                       obj.power_fit.GetParameter(1)),
                       "L")
  probstring = root_utilities.ChisquareString(obj.power_fit)
  legend_prob.AddEntry("",probstring,"")
  legend_prob.Draw()

  raw_input("Displaying graph. Press enter to continue.")
  if save: canvas.SaveAs("plots/power.pdf")

################################################################################

def draw_fourier(obj,save=False):
  x = obj.parent.fourier_x
  y = obj.fourier_y
  n = len(x)
  canvas = TCanvas()
  SetOwnership(canvas,False)
  graph = TGraph(n,x,y)
  SetOwnership(graph,False)
  graph.Draw("AL")
  raw_input("Displaying graph. Press enter to continue.")
  if save: canvas.SaveAs("plots/fourier.pdf")

################################################################################

def draw_height(save=False):

  # Extract data
  height = []
  mean = []
  std = []
  for measurement in anemometer.measurements:
    height.append(measurement.height)
    mean.append(measurement.mean)
    std.append(measurement.std)
  n = len(mean)
  height = np.array(height,dtype="float")
  mean = np.array(mean,dtype="float")
  std = np.array(std,dtype="float")

  # Draw graph
  canvas = TCanvas()
  SetOwnership(canvas,False)
  # graph = TGraph(n,height,mean)
  graph = TGraphErrors(n,height,mean,np.zeros(n),std)
  root_utilities.SetMarker(graph,color=2,style=4,size=2.5)
  graph.SetLineColor(2)
  SetOwnership(graph,False)
  graph.SetTitle(";Height [m];Mean wind speed [m/s]")
  graph.Draw("AP")

  # Fits
  fit_log = TF1("fit_log","[0]*log(([1]*x+[2]))")
  fit_ref = TF1("fit_ref","[0]*(x/[1])^[2]")
  fit_log.SetParameter(1,4.07)
  fit_log.SetParameter(1,0.18)
  fit_log.SetParameter(1,2.22)
  fit_ref.FixParameter(0,mean[0])
  fit_ref.FixParameter(1,height[0])
  root_utilities.SetLine(fit_log,color=8,width=2)
  root_utilities.SetLine(fit_ref,color=9,width=2)
  graph.Fit("fit_log")
  graph.Fit("fit_ref","+")

  # Legend
  legend = TLegend(0.12,0.7,0.5,0.88)
  legend.SetFillStyle(0)
  legend.AddEntry(fit_log,"f_{1}(x) = %.2f*log(%.2fx + %.2f)" % (
                  fit_log.GetParameter(0),
                  fit_log.GetParameter(1),
                  fit_log.GetParameter(2)
                 ),"L")
  legend.AddEntry(fit_ref,"f_{2}(x) = %.2f*(x/%.1f)^{%.2f}" % (
                  fit_ref.GetParameter(0),
                  fit_ref.GetParameter(1),
                  fit_ref.GetParameter(2)
                 ),"L")
  legend.SetLineColor(0)
  legend.Draw()

  # Probability legend
  prob_legend = TLegend(0.4,0.12,0.88,0.4)
  prob_legend.SetFillStyle(0)
  prob_legend.SetLineColor(0)
  fitstring_log = root_utilities.ChisquareString(fit_log)
  fitstring_ref = root_utilities.ChisquareString(fit_ref)
  prob_legend.AddEntry(fit_log,fitstring_log,"L")
  prob_legend.AddEntry(fit_ref,fitstring_ref,"L")
  (_,_,prob_log) = root_utilities.ChisquareStats(fit_log)
  (_,_,prob_ref) = root_utilities.ChisquareStats(fit_ref)
  prob_legend.Draw()

  raw_input("Displaying graph. Press enter to continue.")
  if save: canvas.SaveAs("plots/height.pdf")

################################################################################

folder = "/Users/johannes/Dropbox/Applied Statistics/Project 2/Data/"

anemometer = WindData(folder + "windspeed.data")
sonicmeter = WindData(folder + "sonic.data")

color_reg = []
colors = []
color_index = 1000
reds = np.linspace(0,1,anemometer.n_heights)
for i in range(anemometer.n_heights):
  color_reg.append(TColor(color_index,reds[i],0,1,"",1))
  colors.append(color_index)
  color_index += 1

################################################################################