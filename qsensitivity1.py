import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sst

#Variables
alpha = 0.5
q = np.arange(-1.5,1.51,0.01)
xlst = [1.,5.,10.,20.]
T = 0.02

dKdqlst = []
Klst = []
for x in xlst:
    dKdq = (np.exp(2*q*T)-2*q*T*np.exp(q*T)-1)/(np.exp(q*T)-1)**2
    K = q*(1+np.exp(q*T))/(np.exp(q*T)-1)
    dKdqlst.append(dKdq)
    Klst.append(K)

size = 16
plt.subplot(1,2,1)
plt.plot(q,Klst[0],label='x = '+str(xlst[0])+'m')
#plt.plot(q,Klst[1],'r-',label='x = '+str(xlst[1])+'m')
#plt.plot(q,Klst[2],'k-',label='x = '+str(xlst[2])+'m')
#plt.plot(q,Klst[3],'g-',label='x = '+str(xlst[3])+'m')
plt.xticks(fontsize = size-2)
plt.yticks(fontsize = size-2)
plt.xlabel('q',fontsize=size)
plt.ylabel('K\'(q) [-]',fontsize=size)
plt.legend(loc='upper left')
plt.subplot(1,2,1)


plt.show()

#cov = np.cov(Klst[0],Klst[-1])

#emax = Klst[-1][-1] - Klst[0][-1]


##plt.show()
