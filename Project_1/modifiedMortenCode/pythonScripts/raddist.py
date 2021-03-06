"""
In order to run this program one must have MayaVI installed.
I run it as python raddist.py -toolkit wx
"""
import sys, os, numpy as np

from math import sqrt
from mayavi import mlab
import matplotlib.pyplot as plt
#from pylab import *


def loadCube(positions,N=200, l=-3, u=3):

    bc=np.linspace(l, u, N);
    cube = np.zeros([N,N,N])

    for pos in positions:
        posc=[0,0,0]
        for d in xrange(3):
            for i in xrange(N):
                if not posc[d]:
                    if pos[d] < bc[i]:
                        posc[d] = i
        cube[posc[0],posc[1],posc[2]] += 1

    print "Done filling cube"
    return cube


def radialDistribution(positions, N=100):

    u = positions.max()
    print len(positions)
    r = np.zeros(len(positions))

    for i in xrange(len(positions)):
        r2 = 0.0
        for j in xrange(3):
            r2 += positions[i][j]*positions[i][j]
        r[i] = sqrt(r2);

    print "Done filling vec"
    Weights = np.ones_like(r)/float(len(r))
    n, bins, patches = plt.hist(r, bins=N,normed=1, weights=Weights, color="#006867")
    #plt.hold("on")
    #plt.plot(bins[:-1], n, '--', color="#680000", linewidth=3)
    plt.xlabel("Radial distance")
    plt.ylabel("Probability")
    plt.grid("on")
    plt.show()
    plt.figure()
    plt.plot(bins[:-1], n, '--', color="#680000", linewidth=3)
    plt.grid("on")
    plt.xlabel("Radial distance")
    plt.ylabel("Probability")
    plt.show()






def earthSpherify(data):
    """Creates a spherical representation of the data with a slice to the center"""

    n, m, l = data.shape

    f = lambda i, n: (i*2.0/(n-1) - 1)**2
    f = np.vectorize(f)

    D_i = f(xrange(0, n), n)
    D_j = f(xrange(0, m), m)
    D_k = f(xrange(0, l), l)
    nBins = 100
    rhist = np.linspace(0,1,nBins)
    hist = np.zeros(nBins)
    #Create the sphere
    for i , d_i in enumerate(D_i):
        for j, d_j in enumerate(D_j):
            for k, d_k in enumerate(D_k):
                r = sqrt(d_i + d_j + d_k)
                if r > 1:
                    data[i, j, k] = 0
                else:
                    pos = int(r*nBins)
                    #hist[pos]+=data[i,j,k];
    #Create the slice
    # plt.plot(rhist,hist)
    # plt.show()
    data[n/2:, m/2:, l/2:] = 0

    return data;


def make3dHist(data):
    # n, m, l = data.shape
    # results = np.zeros((n*l,3))
    # for i in range(n):
    #     for j in range(m):
    #         results[i*l:(i+1)*l,j] = data[i,j,:]

    # H, edges = np.histogramdd(results, bins = (100,100,100))
    # print H.shape
    data = earthSpherify(data)
    #mlab.contour3d(data)
    data = data/np.max(data)*2    # Normalize
    mlab.pipeline.volume(mlab.pipeline.scalar_field(data), vmin=0, vmax=1.0)
    mlab.show()


def make1dHist(data):
    Weights = np.ones_like(data)/float(len(data))
    plt.hist(data, bins=np.linspace(0,0.4,33), weights=Weights)
    plt.show()


if __name__ == "__main__":
    path1 = "../datafiles/positionsHOe5.txt"
    path2 = "../datafiles/positionsIe7.txt"
    path3 = "../datafiles/positionsIe6.txt"
    path4 = "../datafiles/positionsHOe6.txt"
    path5 = "../datafiles/positionsIe3.txt"
    path6 = "../datafiles/positionsHOe7.txt"
    positions = np.loadtxt(path6)
    print "Done loading file"


    #make1dHist(data)
    #data = loadCube(positions)
    #make3dHist(data);

    radialDistribution(positions, 100)


    #make2dHist(loadArmaCube(path2), 0)
    #makeRadialHist1(loadArmaCube(path2))
    #makeRadialHist2(loadArmaMat(path6))
    #plt.hold('on')
    #makeRadialHist2(loadArmaMat(path7))
    #makeRadialHist2(loadArmaMat(path8))
    #makeRadialHist2(loadArmaMat(path9))
    # data = loadArmaVec(path3)
    # make1dHist(data);
    # data = loadArmaMat(path4)
    # make2dHist(data);
    # data = loadArmaCube(path5)

    #plt.xlabel('Distance [a.u]')
    #plt.ylabel('Probability')
    #plt.legend(['with correlations', 'without correlations, alpha=10', 'without correlations, alpha=7.81', 'with correlations, alpha=10, beta=0.104'])
    plt.show()
