from math import sqrt
from sys import argv

def std(tau, dt, var):
    return sqrt(1+2*tau/dt)*var

if __name__ == '__main__':
    tau = float(argv[1])
    dt  = float(argv[2])
    var = float(argv[3])
    print std(tau, dt, var)
