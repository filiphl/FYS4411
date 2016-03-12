from numpy import *
from matplotlib import pyplot as plt


N = 5000;

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


blockSizes = range(1,N,1)
std = zeros(len(blockSizes))

index = 0
for blockSize in blockSizes:
	print blockSize
	numberOfBlocks = numberOfSamples/blockSize # Intentional integer division
	averages = zeros(numberOfBlocks)
	newSamples = zeros(numberOfBlocks)

	for i in range(numberOfBlocks):
		lowerBound = i*blockSize
		upperBound = lowerBound + blockSize
		averages[i] = mean(energy[lowerBound:upperBound])

	meanE  = mean(averages)
	meanE2 = mean(averages**2)

	varE = meanE2-meanE**2
	std[index] = sqrt(varE/numberOfBlocks)
	index += 1


plt.plot(blockSizes, std, linewidth=1, color="#1A474A")
plt.xlabel('Block size')
plt.ylabel('Standard deviation')
plt.grid("on")
plt.show()
