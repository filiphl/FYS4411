from numpy import *
from matplotlib import pyplot as plt


N = 400;

infile = open('localenergies.txt', 'r')

energy = []
for line in infile:
	col = line.split()
	energy.append(float(col[0]))

energy = asarray(energy)
numberOfSamples = len(energy)


def standardDeviation(entry):
	e1  = sum(entry)/float(len(entry))
	e2  = sum(entry*entry)/len(entry)
	var = e2 - e1**2
	if var < 0:
		print "Somtn wrong"
	return sqrt(var)

def mean(entry):
	return sum(entry)/len(entry)

averages = []
std = []
for numberOfBlocks in range(1,N):
	print numberOfBlocks
	averages.append([])
	blockSize = numberOfSamples/numberOfBlocks		#intentional integer division
	for block in range(0,numberOfBlocks):
		lowerBound = block*blockSize
		upperBound = lowerBound + blockSize
		averages[numberOfBlocks-1].append(mean(asarray(energy[lowerBound:upperBound])))
	std.append(standardDeviation(asarray(averages[numberOfBlocks-1])))

#for i in xrange(len(sample)):
#	std.append(standardDeviation(asarray(sample[i])))

plt.plot(range(1,N), std, linewidth=3, color="#1A474A")
plt.xlabel('Number of blocks')
plt.ylabel('standard deviation')
plt.show()
