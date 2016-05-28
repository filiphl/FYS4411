from numpy import *
from matplotlib import pyplot as plt
from sys import argv

filename = argv[0]

bs = filename + 'BS.npz'
sd = filename + 'SD.npz'

blockSizes = load(bs)
std        = load(sd)

plt.plot(blockSizes, std, linewidth=1, color="#1A474A")
plt.xlabel('Block size')
plt.ylabel('Standard deviation')
plt.grid('on')
plt.savefig(outFileName, format='pdf')
plt.show()
