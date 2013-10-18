from ROOT import TCanvas, TGraph, TMath, SetOwnership
from numpy import array

# Usage: variable must be defined at caller, as in line = DrawLine(...)
# Examples: line = DrawLine(x1=0,x2=3,y2=2,color=4,width=2)
#           line = DrawLine(x1=-5,x2=5,y=3,style=2,color=2)
#           line = DrawLine(x1=2,x2=4,y1=2,y2=4,color=2,style=2)
def DrawLine(x1=0,x2=0,y1=0,y2=0,x=None,y=None,color=1,width=1,style=1):
  if x != None:
    x1 = x2 = x
  if y != None:
    y1 = y2 = y
  line = TGraph(2,array([x1,x2],dtype="float"),array([y1,y2],dtype="float"))
  SetOwnership(line,False)
  SetLine(line,color=color,width=width,style=style)
  line.Draw("same")
  return line

# Usage: call RemoveAxis on the axis you want to remove.
# Examples: RemoveAxis(my_graph.GetXaxis())
#           RemoveAxis(my_graph.GetYaxis())
def RemoveAxis(axis):
  axis.SetLabelSize(0)
  axis.SetLabelOffset(999)
  axis.SetTickLength(0)
  return axis

# Usage: specify one or more axis limits that should be changed for a plot.
# Examples: SetRanges(my_hist,y2=100)
#           SetRanges(my_graph,x1=0,x2=5,y1=-5,y2=-5)
def SetRanges(plot,x1=None,x2=None,y1=None,y2=None):
  xaxis = plot.GetXaxis()
  yaxis = plot.GetYaxis()
  x1 = x1 if x1 != None else xaxis.GetXmin()
  x2 = x2 if x2 != None else xaxis.GetXmax()
  y1 = y1 if y1 != None else yaxis.GetXmin()
  y2 = y2 if y2 != None else yaxis.GetXmax()
  plot.GetXaxis().SetRangeUser(x1,x2)
  plot.GetYaxis().SetRangeUser(y1,y2)
  return plot

# Usage: specify multiple marker parameters in one function call.
# Examples: SetMarker(my_graph,color=2,style=3,size=1.5)
#           SetMarker(my_graph,color=4,size=3)
#           SetMarker(my_graph,size=2,style=2)
def SetMarker(plot,color=1,style=5,size=2):
  plot.SetMarkerColor(color)
  plot.SetMarkerStyle(style)
  plot.SetMarkerSize(size)
  return plot

# Usage: specify multiple line parameters in one function call.
# Examples: SetLine(my_graph,color=4,style=2)
#           SetLine(my_line,width=2,color=2)
def SetLine(plot,color=1,style=1,width=1):
  plot.SetLineColor(color)
  plot.SetLineStyle(style)
  plot.SetLineWidth(width)
  return plot

# Removes background and border color to make legend stealthy
def StealthyLegend(legend):
  legend.SetFillStyle(0)
  legend.SetLineColor(0)
  return legend

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