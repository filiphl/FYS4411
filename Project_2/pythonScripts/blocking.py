from numpy import *
from matplotlib import pyplot as plt
from sys import argv

filename = str(argv[1])
N 		 = int(argv[2]);
step     = int(argv[3])


outFileName = '../Report/figures/blocking/' + filename[19:-4] + ".pdf"

energy = fromfile(filename)
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


blockSizes = range(step,N,step)
std = zeros(len(blockSizes))

index = 0
for blockSize in blockSizes:
	print '{0}\r'.format(blockSize),
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
plt.grid('on')
#plt.savefig(outFileName, format='pdf')
#save("blockingData/"+filename[19:-4]+"BS",  blockSizes)
#save("blockingData/"+filename[19:-4]+"STD", std)
plt.show()
