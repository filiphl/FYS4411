from numpy import *
from matplotlib import pyplot as plt

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
	return sqrt(var) 


std = []
sample = []
for blockSize in xrange(1,100):
	sample.append([])
	blockSum = 0
	counter = 0
	for e in energy:
		if counter == blockSize:
			sample[blockSize-1].append(float(blockSum)/blockSize)
			blockSum = 0
			counter  = 0
		blockSum += e
		counter += 1

for i in xrange(len(sample)):
	std.append(standardDeviation(asarray(sample[i])))

plt.plot(range(1,100), std)
plt.show()





