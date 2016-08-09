import numpy as np
import matplotlib.pyplot as plt
uz = np.array([-0.791331,   -0.41561585,  0.04095741,  0.56926973,  1.15526113,  1.7798279,  2.41895989,  3.0440846,   3.6225763 ])
wzx = np.array([0.001873434796042591, 0.00081917830353730981, -0.00044509334479308804, -0.0018919171989150573, -0.0034805463311174427, -0.0051568622164867755, -0.0068539012343729003, -0.0084929045570033002, -0.009984773096063838])

uz_bar = np.average(uz)
wzx_bar = np.average(wzx)

diff_uz = uz - uz_bar
diff_wzx = wzx - wzx_bar

cov = sum(diff_uz*diff_wzx)/(uz.shape[0]-1)
print cov

# u_in = np.arange(-15,15,0.01)
# L = 26. #Curve's amplitude
# k = 1.0 #steepness
# x0 = 1. #midway point
#
# # for k in range(1,5):
# #     k = float(k)/2
# u_out = L/(1+np.exp(-k*(u_in-x0))) -L/2
#
# plt.plot(u_in,u_out,label="k="+str(k))
# plt.legend()
# plt.show()