#!/usr/bin/python3
import sys
import math

print('Opening'+sys.argv[1])
f = open(sys.argv[1],'r') # First arg is filename
c = int(sys.argv[2])      # The column holding the values
nSmpl = int(sys.argv[3])  # Number of Samples
bb = int(sys.argv[4])     # Largest bin bin_size

valAvg = 0
valLst = []
# Read first two lines - header + blank line
f.readline()
f.readline()

eql = 0
prl = 0
for i in range(1,nSmpl+1):
	fields = f.readline().split()
	prl = prl+1
	valAvg = valAvg + float(fields[c])
	valLst.append(float(fields[c]))
print('Num of samples: '+str(prl))
valAvg = valAvg / prl
print('Val_smpl_avg: {0:.5f}'.format(valAvg))

#Start binning
print('Bin size, StdDev')
val_smpl_var = 0.0

for bin_size in range(1,bb+1):
	# We will loop only through whole number bb-tuples
	# discarding the rest
	val_b = []
	d_point = iter(valLst)
	bin_num = prl//bin_size
	for bin_no in range(bin_num):
		val_bin_avg = 0.0
		for d in range(bin_size):
			val_bin_avg += next(d_point)
		val_bin_avg = val_bin_avg/bin_size
		val_b.append(val_bin_avg)
	val_tot_var = 0.0
	for bin_no in range(bin_num):
		val_tot_var += (val_b[bin_no]-valAvg)**2.0
	val_tot_var = val_tot_var/bin_num
	if(bin_size == 1):
		# Computes and sets the variance for bin_size = 1
		val_smpl_var = val_tot_var
	val_mean_var = math.sqrt(val_tot_var/bin_num)
	print('{0:05d} {1:.5f} {2:.5f} {3:.5f}'.format( \
		bin_size, val_mean_var, math.sqrt(bin_size), \
		bin_size*val_tot_var/val_smpl_var))