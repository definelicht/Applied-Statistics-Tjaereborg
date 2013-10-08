from numpy import ndarray, zeros
import csv

# Examples: data = ReadCsv("mydata.txt",delimiter="\t")
#           data = ReadCsv("nicedata.csv",skip=1)
def ReadCsv(filename,delimiter=",",skip=0):
  datafile = open(filename,"rU")
  reader = csv.reader(datafile,delimiter=delimiter)
  col_count = len(reader.next())
  row_count = sum(1 for row in reader) + 1
  datafile.seek(0)
  for i in range(skip):
    reader.next()
    row_count -= 1
  output = zeros((row_count,col_count),dtype="float")
  i = 0
  for row in reader:
    output[i] = row
    i += 1
  datafile.close()
  return output

def Header(filename,delimiter=","):
  datafile = open(filename,"rU")
  reader = csv.reader(datafile,delimiter=delimiter)
  header = reader.next()
  datafile.close()
  return header