from numpy import *
from matplotlib import pyplot as plt

infile = open("../datafiles/positions.txt")

positions = []
for line in infile:
    col = line.split()
    positions.append([float(col[0]), float(col[1]), float(col[2])])

maximum = 0;
minimum = 0;
for p in positions:
    if max(p)>maximum:
        maximum = max(p)
    if min(p) < minimum:
        minimum = min(p)

print "max: ", maximum
print "min: ", minimum

N = 100
grid = linspace(-4,4,N)
cells = zeros([N,N,N])

for pos in positions:
    x = 0; y=0; z=0;
    for i in xrange(N):
        if pos[0]<grid[i]:
            x=i
        if pos[1]<grid[i]:
            y=i
        if pos[2]<grid[i]:
            z=i
    cells[x,y,z] += 1
