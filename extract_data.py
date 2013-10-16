import sys, numpy as np, read_csv

# Usage: python extract_data.py folder=/jesus/applied_statistics/data/

folder = "/users/johannes/Dropbox/Applied Statistics/Project 2/Data/"
for arg in sys.argv:
  if arg[0:5] == "folder=":
    folder = arg[5:]
filename = folder + "all.data"

data = read_csv.ReadCsv(filename)
lines = len(data)
windspeed_file = open(folder + "windspeed.data","w")
windspeed_file.write("t,h=17m,h=28.5m,h=41m,h=57m,h=77m,h=90m")
sonic_file = open(folder + "sonic.data","w")
sonic_file.write("t,h=17m,h=77m,h=90m")
for i in range(lines):
  windspeed_file.write("\n%f,%f,%f,%f,%f,%f,%f" % (data[i,173],data[i,117],data[i,116],data[i,115],data[i,114],data[i,113],data[i,112]))
  sonic_file.write("\n%f,%f,%f,%f" % (data[i,173],data[i,109],data[i,101],data[i,93]))
windspeed_file.close()
sonic_file.close()