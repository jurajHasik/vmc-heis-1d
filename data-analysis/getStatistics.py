import math, sys

dataFile = sys.argv[1] # File containing data
nSmpl    = int(sys.argv[2]) # number of samples in dataFile 

f = open(dataFile, 'r')
# Read first two lines - header + blank line
f.readline()
f.readline()

eloc = []
for i in range(1,nSmpl):
	line = f.readline().split()
	# 2nd column (1st index in line) holds e_loc
	eloc.append(float(line[1]))

eAvg = 0.0
for e in eloc:
	eAvg += e
eAvg = eAvg/nSmpl

eVar = 0.0
for e in eloc:
	eVar += (e-eAvg)**2.0
eVar = math.sqrt(eVar/nSmpl)

print("AVG(e_loc)= {0:10.5f} STDEV(e_loc)= {1:10.5f}"\
	.format(eAvg, eVar))