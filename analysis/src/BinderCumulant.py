# binder.py
# 
# A script to compute the binder cumulant
#

import sys
import math

args = sys.argv

if (len(args) < 8):
    print "Usage: binder.py [time_col] [value_col] [tstart] [tend] [freq] [output_file] [data_files]"
    sys.exit(1)

args.pop(0) # Ignore self
time_col = int(args.pop(0))
value_col = int(args.pop(0))
tstart = int(args.pop(0))
tend = int(args.pop(0))
freq = int(args.pop(0))
output_file = args.pop(0)

max_col = 0
if (time_col > value_col):
    max_col = time_col
else:
    max_col = value_col

# open data files
files = [open(i, "r") for i in args]

count = 0

avg2 = 0.0
avg4 = 0.0

for f in files:
    for line in f:
        if (not line.startswith("#")):
            data = line.strip().split()
            if (data == [] or len(data) < (max_col+1)):
                continue
            time = int(data[time_col])
            if (time < tstart or time > tend or time % freq != 0):
                continue
            value = float(data[value_col])
            avg2 += value ** 2
            avg4 += value ** 4
            count += 1
    f.close()

avg2 /= float(count)
avg4 /= float(count)

bin_cum = 1 - (avg4 / (3 * avg2 * avg2))

writer = open(output_file, "w")
output = "%.5f\n" % bin_cum
writer.write(output)
