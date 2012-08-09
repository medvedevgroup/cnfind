import matplotlib
matplotlib.use('Agg')
import sys
import numpy as np
#from pylab import *
import matplotlib.pyplot as plt
import math

print "hello\n"
width = 0.5 
linetypes = ('k-v', 'r-s', 'g-^', 'b-o', 'c-D', 'm-*', 'k-', 'k-', 'k-', 'k-', 'k-', 'k-')
datafile = sys.argv[1]
legendloc = 4 #bottom right #4 bottom right, 7 center right, 1 upper right


plt.subplots_adjust(left=.15, bottom=.2, right=0.85, top=.9, wspace=None, hspace=None)

f = open(datafile, "r")
x = []
y  = []
count = 0


for line in f:
	fields = line.split()
	if len(fields) == 2 : 
		x.append(fields[0])
		y.append(fields[1])
		count +=1

plt.plot(x,y, markersize=10, linewidth = width )
#plt.plot(x,y,linetypes[0], markersize=10, linewidth = width )


plt.xlabel('Position', fontsize=30 )
plt.ylabel('Data', fontsize=30)
tickfontsize = 20
plt.yticks(size=tickfontsize)
plt.xticks(size=tickfontsize)

leg = plt.legend(loc=legendloc) #4 bottom right, 7 center right, 1 upper right

plt.savefig(datafile + ".ps")
