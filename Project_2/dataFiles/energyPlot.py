from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from numpy import *

data = loadtxt('energyPlotfile.txt')

alpha  = data[:,0]
beta   = data[:,1]
energy = data[:,2]


print alpha
print beta
print energy

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(alpha, beta, energy,  cmap='hot', vmin=min(energy), vmax=max(energy))
plt.xlabel("alpha")
plt.ylabel("beta")
ax.set_zlabel("energy")
plt.show()
