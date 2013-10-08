from __future__ import division
import sys, numpy as np, read_csv
from ROOT import TCanvas, TGraph, TMultiGraph, TLegend, TF1

def WeightedMean(values,errors):
  n = len(values)
  average = 0
  norm = 0
  for i in range(n):
    average += values[i] / errors[i]**2
    norm += errors[i]**(-2)
  average /= norm
  error = 0
  for i in range(n):
    error += errors[i]**(-2)
  error = error**(-0.5)
  return (average,error)

folder = "/users/johannes/Dropbox/Applied Statistics/Project 2/Data/"
hold = False
for arg in sys.argv:
  if arg[0:7] == "folder=":
    folder = arg[7:]
  if arg == "hold":
    hold = True
filename = folder + "windspeed.data"

header = read_csv.Header(filename)
height = []
for h in header:
  try:
    height.append(float(h[2:len(h)-1]))
  except:
    continue
height = np.array(height,dtype="float")
n_heights = len(height)

windspeed = read_csv.ReadCsv(filename,skip=True)
time = np.array(windspeed[:,0])
n = len(time)
windspeed = windspeed[:,1:]
t_lim = [np.min(time),np.max(time)]
wind_lim = [np.min(windspeed),np.max(windspeed)]

multigraph = TMultiGraph()
graphs = []
gauss = []
legend = TLegend(0.8,0.7,0.95,0.95)
mean = np.zeros(n_heights,dtype="float")
std = np.zeros(n_heights,dtype="float")
for i in range(n_heights):
  # Plot raw data
  wind = np.array(windspeed[:,i])
  graph = TGraph(n,time,wind)
  graph.GetYaxis().SetRangeUser(wind_lim[0],wind_lim[1])
  graph.SetLineColor(i+1)
  legend.AddEntry(graph,"h = %.1fm" % height[i-1],"L")
  graphs.append(graph)
  multigraph.Add(graph)
  # Calculate normal distribution
  mean[i] = np.mean(wind)
  std[i] = np.std(wind)
  f = TF1("gaus%i"%i,"gaus",wind_lim[0],wind_lim[1])
  f.SetParameters(height[i],mean[i],std[i])
  f.SetLineColor(i+1)
  gauss.append(f)

# Aggregate analysis
mean_fit1 = TF1("mean_fit1","[0]+[1]*x+[2]*x^2")
mean_fit2 = TF1("mean_fit2","[0]*exp([1]*x)")
mean_fit3 = TF1("mean_fit3","[0]*pow(x,[1])")
mean_fit1.SetLineColor(6)
mean_fit2.SetLineColor(8)
mean_fit3.SetLineColor(9)
graph_mean = TGraph(n_heights,height,mean)
graph_mean.SetMarkerStyle(2)
graph_mean.SetMarkerColor(4)
graph_mean.SetMarkerSize(2.5)
graph_mean.SetTitle(";Height [m];Mean wind speed [m/s^{2}]")
canvas_mean = TCanvas()
graph_mean.Draw("AP")
graph_mean.Fit("mean_fit1")
graph_mean.Fit("mean_fit2","+")
graph_mean.Fit("mean_fit3","+")
mean_legend = TLegend(0.14,0.7,0.4,0.85)
mean_legend.AddEntry(mean_fit1,"%.2fx^{2} + %.2fx - %.2f" % (mean_fit1.GetParameter(0),mean_fit1.GetParameter(1),abs(mean_fit1.GetParameter(2))),"L")
mean_legend.AddEntry(mean_fit2,"Exponential","L")
mean_legend.AddEntry(mean_fit3,"Power","L")
mean_legend.SetFillColor(0)
mean_legend.SetLineColor(0)
mean_legend.Draw()

# Draw graphs
canvas_multigraph = TCanvas()
multigraph.Draw("AL")
multigraph.SetTitle(";Time [s];Wind speed [m/s^{2}]")
legend.SetFillColor(0)
legend.Draw()

canvas_multigraph.SaveAs("plots/rawdata.pdf")

if hold: raw_input("Holding...")