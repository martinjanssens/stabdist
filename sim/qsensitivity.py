import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import scipy.stats as sst

#Variables
alpha = 0.5
qlst = [-1.5,-0.1,0.1,1.5]
x = np.arange(0,10.1,0.1)
T = 0.02

dKlst = []
Klst = []
for q in qlst:
    dKdq = (np.exp(2*q*T)-2*q*T*np.exp(q*T)-1)/(np.exp(q*T)-1)**2*x
    K = q*(1+np.exp(q*T))/(np.exp(q*T)-1)*x
    dKlst.append(dKdq)
    Klst.append(K)

qrange = np.arange(-1.5,1.51,0.01)
Kq = qrange*(1+np.exp(qrange*T))/(np.exp(qrange*T)-1)

size = 22
plt.subplot(1,2,1)
plt.plot(qrange,Kq)#,label='x = '+str(xlst[0])+'m')
# plt.xticks(fontsize = size-2)
# plt.yticks(fontsize = size-2)
plt.xlabel(r'q [1/s]')#,fontsize=size)
plt.ylabel('K/x [1/m]')#,fontsize=size)
ax = plt.gca()
ax.set_xlim([-1.5,1.5])
ax.get_yaxis().get_major_formatter().set_useOffset(False)
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
# plt.legend(loc='upper left')
plt.subplot(1,2,2)
plt.plot(x,Klst[0],'b-',label='q = '+str(qlst[0]))
plt.plot(x,Klst[1],'k-',label='q = '+str(qlst[1]))
plt.plot(x,Klst[2],'r-',label='q = '+str(qlst[2]))
plt.plot(x,Klst[3],'g-',label='q = '+str(qlst[3]))
# plt.xticks(fontsize = size-2)
# plt.yticks(fontsize = size-2)
plt.xlabel('x [m]')#,fontsize=size)
plt.ylabel('K [-]')#,fontsize=size)
plt.legend(loc='upper left', fontsize=20)
matplotlib.rcParams.update({'font.size': 20})  # increase font size on axes (edited)
plt.savefig('output\qsens.eps', bbox_inches = 'tight')
plt.show()

cov = np.cov(Klst[0],Klst[-1])

emax = Klst[-1][-1] - Klst[0][-1]


##plt.show()
