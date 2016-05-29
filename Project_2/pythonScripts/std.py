from math import sqrt
from sys import argv

def std(tau, var):
    return sqrt((1+2*tau)*var)

if __name__ == '__main__':
    tau = float(argv[1])
    var = float(argv[2])
    print std(tau, var)
