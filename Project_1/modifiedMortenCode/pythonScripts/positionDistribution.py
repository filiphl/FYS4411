from numpy import *
from matplotlib import pyplot as plt
import sys
infile = open("../datafiles/positions.txt")
"""
positions = []
for line in infile:
    col = line.split()
    positions.append([float(col[0]), float(col[1]), float(col[2])])
"""
positions = loadtxt("../datafiles/positions.txt")
numberOfSteps = len(positions)

N = 1000
grid = linspace(-2,2,N)
cells = zeros([N,N])

counter = 0;
for pos in positions:
    sys.stdout.write("  %d%% complete\r"%(100.0*counter/numberOfSteps))
    sys.stdout.flush()
    x = 0; y=0; z=0;
    for i in xrange(N):
        if not x:
            if pos[0]<grid[i]:
                x=i
        if not y:
            if pos[1]<grid[i]:
                y=i
        if (x and y):
            break
    counter += 1
    cells[x,y] += 1


plt.imshow(cells)
plt.show()
