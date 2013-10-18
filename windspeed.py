from __future__ import division
import math, sys, numpy as np, read_csv, root_utilities, os
from ROOT import TCanvas, TGraph, TGraphErrors, TH1F, TMultiGraph, TLegend, TF1, TMath, TPad, TColor

def GaussFromHist(data,nbins=None):
  if nbins == None: nbins = int(len(data)/20)
  hist = TH1F(os.urandom(8),"",nbins,min(data),max(data))
  for i in data: hist.Fill(i)
  fitname = os.urandom(8)
  fit = TF1(fitname,"gaus")
  hist.Fit(fitname)
  return fit

""" Parse command line arguments """
plot_prefix = "plots/"
folder = "/Users/johannes/Dropbox/Applied Statistics/Project 2/Data/"
hold = False
single = -1
sonic = False
for arg in sys.argv:
  if arg[:7] == "folder=":
    folder = arg[7:]
  if arg == "hold":
    hold = True
  if arg[:7] == "single=":
    single = int(arg[7:])
  if arg == "save":
    plot_prefix = "/Users/johannes/Dropbox/Applied Statistics/Project 2/Figurer/"
  if arg == "sonic":
    sonic = True
if sonic:
  filename = folder + "sonic.data"
  plot_prefix += "sonic_"
else:
  filename = folder + "windspeed.data"


""" Read header and parse heights """
header = read_csv.Header(filename)
height = []
for h in header:
  try:
    height.append(float(h[2:len(h)-1]))
  except:
    continue
height = np.array(height,dtype="float")
n_heights = len(height)


""" Read data and set general variables """
windspeed = read_csv.ReadCsv(filename,skip=True)
time = np.array(windspeed[:,0])
n = len(time)
windspeed = windspeed[:,1:]
t_lim = [np.min(time),np.max(time)]
wind_lim = [np.min(windspeed),np.max(windspeed)]


""" Setup colors """
color_reg = []
colors = []
reds = np.linspace(0,1,n_heights)
for i in range(n_heights):
  color_reg.append(TColor(i+1000,reds[i],0,1,"",1))
  colors.append(i+1000)


""" Treat data """
rawdata_multigraph = TMultiGraph()
fourier_multigraph = TMultiGraph()
multipower_multigraph = TMultiGraph()
graphs = []
fourier_graphs = []
# gauss = []
rawdata_legend = TLegend(0.8,0.7,0.95,0.95)
fourier_legend = TLegend(0.7,0.55,0.88,0.88)
mean = np.zeros(n_heights,dtype="float")
std = np.zeros(n_heights,dtype="float")
next_power = 2**(int(math.ceil(math.log(n, 2))))
fourier_x = 0.5*35.0*np.linspace(0,1,next_power>>1)
power_y = 0
for i in range(n_heights):
  # Plot raw data
  wind = np.array(windspeed[:,i])
  """  Fourier junk """
  wind_fourier = np.fft.fft(wind,next_power>>1) / n
  fourier_y = 2.0*np.absolute(wind_fourier[:next_power>>1])
  if i == n_heights-1: power_y = fourier_y
  fourier_graph = TGraph(len(fourier_x),fourier_x,fourier_y)
  fourier_graph.SetLineColor(colors[i])
  fourier_graphs.append(fourier_graph)
  if single == -1 or single == i:
    fourier_multigraph.Add(fourier_graph)
    multipower_multigraph.Add(fourier_graph)
    fourier_legend.AddEntry(fourier_graph,"h = %.1fm" % height[i-1],"L")
  """ /Fourier junk """
  graph = TGraph(n,time,wind)
  graph.GetYaxis().SetRangeUser(wind_lim[0],wind_lim[1]+5)
  graph.SetLineColor(i+1)
  rawdata_legend.AddEntry(graph,"h = %.1fm" % height[i-1],"L")
  graphs.append(graph)
  rawdata_multigraph.Add(graph)
  # Calculate normal distribution
  gausfit = GaussFromHist(wind)
  mean[i] = gausfit.GetParameter(1)
  std[i] = gausfit.GetParameter(2)
  # f = TF1("gaus%i"%i,"gaus",wind_lim[0],wind_lim[1])
  # f.SetParameters(height[i],mean[i],std[i])
  # f.SetLineColor(colors[i])
  # gauss.append(f)


""" Aggregate """
mean_fit1 = TF1("mean_fit1","[0]+[1]*x+[2]*x^2")
# mean_fit2 = TF1("mean_fit2","[0]*exp([1]*x)")
mean_fit3 = TF1("mean_fit3","[0]*pow(x,[1])")
mean_fit1.SetLineColor(6)
# mean_fit2.SetLineColor(8)
mean_fit3.SetLineColor(9)
mean_graph = TGraphErrors(n_heights,height,mean,np.zeros(n_heights),std)
mean_graph.SetMarkerStyle(3)
mean_graph.SetMarkerSize(3)
mean_graph.SetTitle(";Height [m];Mean wind speed [m/s^{2}]")
mean_canvas = TCanvas()
mean_graph.Draw("AP")
mean_graph.Fit("mean_fit1")
# mean_graph.Fit("mean_fit2","+")
mean_graph.Fit("mean_fit3","+")
mean_legend = TLegend(0.11,0.70,0.58,0.89)
fits_legend = TLegend(0.55,0.11,0.89,0.30)
chisquares = [mean_fit1.GetChisquare(),mean_fit3.GetChisquare()]
ndf = [mean_fit1.GetNDF(),mean_fit3.GetNDF()]
mean_legend.AddEntry(mean_fit1,root_utilities.ChisquareString(mean_fit1),"L")
# mean_legend.AddEntry(mean_fit2,root_utilities.ChisquareString(mean_fit2),"L")
mean_legend.AddEntry(mean_fit3,root_utilities.ChisquareString(mean_fit3),"L")
pars1 = mean_fit1.GetParameters()
# pars2 = mean_fit2.GetParameters()
pars3 = mean_fit3.GetParameters()
fits_legend.AddEntry(mean_fit1,"%.2fx^{2} + %.2fx - %.2f" % (pars1[0],pars1[1],abs(pars1[2])),"L")
# fits_legend.AddEntry(mean_fit2,"%.2fe^{%.2fx}" % (pars2[0],pars2[1]),"L")
fits_legend.AddEntry(mean_fit3,"%.2fx^{%.2f}" % (pars3[0],pars3[1]),"L")
mean_legend.SetFillColor(0)
mean_legend.SetLineColor(0)
fits_legend.SetFillStyle(0)
fits_legend.SetLineColor(0)
mean_legend.Draw()
fits_legend.Draw()


""" Raw data """
rawdata_canvas = TCanvas()
rawdata_multigraph.Draw("AL")
rawdata_multigraph.SetTitle(";Time [s];Wind speed [m/s^{2}]")
rawdata_legend.SetFillColor(0)
rawdata_legend.Draw()
rawdata_line = []
for i in range(n_heights):
  rawdata_line.append(root_utilities.DrawLine(x1=t_lim[0],x2=t_lim[1],y=mean[i],color=i+1,style=2,width=3))
  rawdata_line[i].Draw("same")


""" Fourier transform """
fourier_canvas = TCanvas()
fourier_multigraph.Draw("AL")
fourier_multigraph.SetTitle(";f [Hz];|Y(f)|")
fourier_multigraph.GetXaxis().SetLimits(0,0.5)
fourier_multigraph.SetMinimum(0)
fourier_multigraph.SetMaximum(0.3)
fourier_legend.SetFillStyle(0)
fourier_legend.SetLineColor(0)
fourier_legend.Draw()
fourier_canvas.Update()


""" Power fit """
power_canvas = TCanvas()
power_canvas.SetLogx()
power_canvas.SetLogy()
power_cutoff = [2e-1,2e0]
power_fit = TF1("powerfit","[0]*x^[1]",power_cutoff[0],power_cutoff[1])
power_fit.SetParameter(0,0.01)
power_fit.SetParameter(1,-1.1)
root_utilities.SetLine(power_fit,style=1,color=2,width=4)
power_x1 = []
power_x2 = []
power_x3 = []
power_y1 = []
power_y2 = []
power_y3 = []
n_fourier = len(power_y)
for i in range(n_fourier):
  x = fourier_x[i]
  y = power_y[i]
  if fourier_x[i] < power_cutoff[0]:
    power_x1.append(x)
    power_y1.append(y)
  elif fourier_x[i] > power_cutoff[1]:
    power_x3.append(x)
    power_y3.append(y)
  else:
    power_x2.append(x)
    power_y2.append(y)
power_x1 = np.array(power_x1,dtype="float")
power_x2 = np.array(power_x2,dtype="float")
power_x3 = np.array(power_x3,dtype="float")
power_y1 = np.array(power_y1,dtype="float")
power_y2 = np.array(power_y2,dtype="float")
power_y3 = np.array(power_y3,dtype="float")
power_graph = TMultiGraph()
power_graph1 = TGraph(len(power_x1),power_x1,power_y1)
power_graph1.SetLineColor(4)
power_graph2 = TGraph(len(power_x2),power_x2,power_y2)
power_graph2.SetLineColor(7)
power_graph3 = TGraph(len(power_x3),power_x3,power_y3)
power_graph3.SetLineColor(4)
power_graph.Add(power_graph1)
power_graph.Add(power_graph2)
power_graph.Add(power_graph3)
power_graph.Draw("AL")
power_graph.SetTitle(";f [Hz];|Y(f)|")
power_graph.GetXaxis().SetLimits(1e-2,1e1)
power_graph.GetYaxis().SetRangeUser(1e-6,1e0)
power_graph2.Fit("powerfit","r")
power_legend = TLegend(0.5,0.7-0.5,0.88,0.88-0.5)
power_legend.SetFillStyle(0)
power_legend.SetLineColor(0)
power_legend.AddEntry(power_graph2,"Fitted data","L")
power_legend.AddEntry(power_fit,"f(x) = %.2ex^{%.2f}"%(power_fit.GetParameter(0),power_fit.GetParameter(1)),"L")
power_chisquare = root_utilities.ChisquareStats(power_fit)
power_legend.AddEntry("","#chi^{2} = %.2e, P = %.2e"%(power_chisquare[0],power_chisquare[2]),"")
power_legend.Draw()


""" Power multigraph """
multipower_canvas = TCanvas()
multipower_canvas.SetLogx()
multipower_canvas.SetLogy()
multipower_multigraph.Draw("AL")
multipower_multigraph.SetTitle(";Frequency [Hz];Amplitude [arbitrary]")
multipower_multigraph.GetXaxis().SetLimits(1e-3,1e1)
multipower_multigraph.SetMinimum(1e-6)
fourier_legend.Draw()
multipower_canvas.Update()


""" Save plots """
rawdata_canvas.SaveAs(plot_prefix + "rawdata.pdf")
fourier_canvas.SaveAs(plot_prefix + "fourier.pdf")
power_canvas.SaveAs(plot_prefix + "power.pdf")
multipower_canvas.SaveAs(plot_prefix + "multipower.pdf")
mean_canvas.SaveAs(plot_prefix + "mean.pdf")

if hold: raw_input("Holding...")