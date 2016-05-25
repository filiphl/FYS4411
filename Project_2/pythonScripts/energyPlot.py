from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from numpy import *


y      = linspace(0.2, 3, 11)
x      = linspace(0.1, 1, 11)
energy = matrix(loadtxt('../dataFiles/energyPlotfileN6.txt'))
alpha, beta = meshgrid(x,y)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(alpha, beta, energy, cstride=1, rstride=1, cmap='coolwarm')
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(alpha, beta, energy,  cmap='hot', vmin=min(energy), vmax=max(energy))
plt.xlabel(r"$\beta$", fontsize=18)
plt.ylabel(r"$\alpha$", fontsize=18)
plt.xticks(linspace(0.1,1.1,6))
ax.set_zlabel(r"$\langle E_L \rangle $", fontsize=14, rotation=170)
ax.view_init(15, 170)
#plt.savefig("energyPlot15x15N2.pdf", transparent=True)
plt.show()
