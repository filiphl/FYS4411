import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def vector(path):
  out = np.zeros(5)
  for i,name in enumerate(['1','2','4','8','16']):
    temp = np.loadtxt(path+'/'+name+'.out')   
    out[i] = np.mean(temp)/1e6
  return out

def printout(a):
  for j in a:
    for i in j:
      print(i)
    print('\n')


unvec = vector('unvec')
vecO3 = vector('vecO3')
vecO2 = vector('vecO2')
ratio2 = unvec/vecO2
ratio3 = unvec/vecO3
parallel = unvec[0]/np.array([unvec[0],unvec[1],unvec[2],unvec[3],unvec[4]])


def plot(show=1):
  fig = plt.figure(figsize=(11,7))
  gs = gridspec.GridSpec(2,2)
  gs.update(hspace=0,wspace=0.001)
  ax1 = plt.subplot(gs[0,0])
  ax2 = plt.subplot(gs[1,0])
  ax3 = plt.subplot(gs[:,-1])
  x = [1,2,4,8,16]
  ax1.plot(x,unvec,marker='o',label='Non vectorized')
  ax1.plot(x,vecO2,marker='o',label='Optimization O2')
  ax1.plot(x,vecO3,marker='o',label='Optimization O3')
  ax2.plot(x,ratio2,marker='^',label='Speed up with O2')
  ax2.plot(x,ratio3,marker='^',label='Speed up with O3')
  ax3.plot(x,parallel,marker='o',label='Parallel speed up')
  ax1.set_ylabel('$T$ [s]')
  ax1.set_xlabel('$p$')
  ax2.set_xlabel('$p$')
  ax3.set_xlabel('$p$')
  ax2.set_ylabel('$T_{O0}/T_{OX}$')
  ax3.set_ylabel('$T_s/T_p$')
  ax1.legend()
  ax2.legend()
  ax3.legend()
  #plt.setp(ax1.get_xticklabels(),visible=False)
  gs.tight_layout(fig)
  if show:
    plt.show()
  else:
    plt.savefig('figure1.pdf')


#print([unvec,vecO2,vecO3,ratio2,ratio3])
plot(0)

