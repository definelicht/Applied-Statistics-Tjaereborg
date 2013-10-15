from __future__ import division
import math, sys, numpy as np, read_csv
from ROOT import TCanvas, TGraph, TGraphErrors, TMultiGraph, TLegend, TF1, TMath, TPad

def ChisquareStats(fit):
  chisquare = fit.GetChisquare()
  ndf = fit.GetNDF()
  probability = TMath.Prob(chisquare,ndf)
  return (chisquare,ndf,probability)

def ChisquareString(fit):
  return "#chi^{2} = %.2e, NDF = %i, p = %.2e" % ChisquareStats(fit)

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

rawdata_multigraph = TMultiGraph()
fourier_multigraph = TMultiGraph()
graphs = []
fourier_graphs = []
gauss = []
legend = TLegend(0.8,0.7,0.95,0.95)
mean = np.zeros(n_heights,dtype="float")
std = np.zeros(n_heights,dtype="float")
for i in range(n_heights):
  # Plot raw data
  wind = np.array(windspeed[:,i])
  wind_fourier = np.fft.fft(wind)
  next_power = 2**(int(math.ceil(math.log(n, 2))))
  graph = TGraph(n,time,wind)
  graph.GetYaxis().SetRangeUser(wind_lim[0],wind_lim[1]+5)
  graph.SetLineColor(i+1)
  legend.AddEntry(graph,"h = %.1fm" % height[i-1],"L")
  graphs.append(graph)
  rawdata_multigraph.Add(graph)
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
graph_mean = TGraphErrors(n_heights,height,mean,np.zeros(n_heights),std)
graph_mean.SetMarkerStyle(3)
graph_mean.SetMarkerSize(3)
graph_mean.SetTitle(";Height [m];Mean wind speed [m/s^{2}]")
canvas_mean = TCanvas()
graph_mean.Draw("AP")
graph_mean.Fit("mean_fit1")
graph_mean.Fit("mean_fit2","+")
graph_mean.Fit("mean_fit3","+")
mean_legend = TLegend(0.11,0.7,0.58,0.89)
fits_legend = TLegend(0.55,0.11,0.89,0.30)
chisquares = [mean_fit1.GetChisquare(),mean_fit2.GetChisquare(),mean_fit3.GetChisquare()]
ndf = [mean_fit1.GetNDF(),mean_fit2.GetNDF(),mean_fit3.GetNDF()]
mean_legend.AddEntry(mean_fit1,ChisquareString(mean_fit1),"L")
mean_legend.AddEntry(mean_fit2,ChisquareString(mean_fit2),"L")
mean_legend.AddEntry(mean_fit3,ChisquareString(mean_fit3),"L")
pars1 = mean_fit1.GetParameters()
pars2 = mean_fit2.GetParameters()
pars3 = mean_fit3.GetParameters()
fits_legend.AddEntry(mean_fit1,"%.2fx^{2} + %.2fx - %.2f" % (pars1[0],pars1[1],abs(pars1[2])),"L")
fits_legend.AddEntry(mean_fit2,"%.2fe^{%.2fx}" % (pars2[0],pars2[1]),"L")
fits_legend.AddEntry(mean_fit3,"%.2fx^{%.2f}" % (pars3[0],pars3[1]),"L")
mean_legend.SetFillStyle(0)
mean_legend.SetLineColor(0)
fits_legend.SetFillStyle(0)
fits_legend.SetLineColor(0)
mean_legend.Draw()
fits_legend.Draw()

# Draw graphs
canvas_multigraph = TCanvas()
rawdata_multigraph.Draw("AL")
rawdata_multigraph.SetTitle(";Time [s];Wind speed [m/s^{2}]")
legend.SetFillColor(0)
legend.Draw()
canvas_multigraph.SaveAs("plots/rawdata.pdf")

# Fourier transform


if hold: raw_input("Holding...")