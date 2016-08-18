##############################
# Experimental data analysis #
##############################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import scipy.stats as sst

# Which file to input
bat = str(3)
day = str(2)
run = str(14)
data_name = "bat"+str(bat)+"_day"+str(day)+"_run"+str(run)
data_file = "data/"+data_name+".data"

# Open the .data file with all the messages
raw_data = open(data_file, "r")
raw_data_lines = raw_data.readlines()

# Extract only the FAKE_OPTICAL_FLOW message
data_lines = []
for i in xrange(len(raw_data_lines)):
    # Split the lines
    raw_data_lines[i] = raw_data_lines[i].split()
    # Only append lines with the FAKE_OPTICAL_FLOW label
    if raw_data_lines[i][2] == "FAKE_OPTICAL_FLOW":
        # Convert the data from strings to floats, expect for the label
        for j in range(len(raw_data_lines[i])):
            if j != 2:
                raw_data_lines[i][j] = float(raw_data_lines[i][j])
        # Delete the message label so the rest can be converted to floats
        del raw_data_lines[i][2]
        # Append
        data_lines.append(raw_data_lines[i])

# Convert to array for further processing
data_arr = np.asarray(data_lines)

# Cols: 0 - Time, 1 - sender_id, 2 - optical_flow, 3 - optical_flow_vision_dt (not used), 4 - normalized_thrust,
# 5 - cov_opticalflow, 6 - pstate, 7 - pused, 8 - agl, 9 - x_pos

t_start = 91. #s
found_start = False
t_stop = 116.5 #s
found_stop = False
for t in range(len(data_arr[:,0])):
    if data_arr[:,0][t] > t_start and found_start == False:
        t_start_index = t
        found_start = True
    elif data_arr[:,0][t] > t_stop and found_stop == False:
        t_stop_index = t
        found_stop = True
    elif t == len(data_arr[:,0])-1 and not found_stop:
        t_stop_index = t

# Find average K_zx and std in the (manually identified) instability time range, assuming

K_avg = np.mean(data_arr[t_start_index:t_stop_index,7])
K_std = np.std(data_arr[t_start_index:t_stop_index,7])

print r'$K_{avg} = $', K_avg
print r'$\sigma = $', K_std

# Processed data of day 1 runs:
# day_1_x = [10.4765, 10.4765, 10.4765, 8.762, 9.941, 9.910]
# day_1_K = [134.572, 187.116, 86.478, 48.926, 118.125, 99.826]

# Processed data of day 2 runs:
day_2_x = np.array([5.373,5.367,5.374,3.762,3.762,3.762,6.95,6.95,6.94,4.66,4.665,4.68,6.12,6.11,6.115])
day_2_K_avg = np.array([258403.883523,252827.188368,260810.735937,125325.322004,123363.,
                        135741.,437113.,463948.373832,457853.311719,197835.774306,202602.173032,
                        205199.515234,366969.585258,363208.785301,366392.011814])
day_2_K_std = np.array([4462.44036676,4281.62543838,3589.50378622,7759.17875367,4490.,
                        6154.,  20225., 24114.9739642,29990.0241012,8157.51820668,10918.8854273,
                        11817.3242842,14458.4202971,11951.0054958,12935.6802922])

slope, intercept, r, pvalue, std = sst.linregress(day_2_x, day_2_K_avg)
x_range = np.arange(min(day_2_x), max(day_2_x), 0.01)
K_range = slope * x_range + intercept

print "p-value = ", pvalue
print "r^2 = ", r*r

# Plotting
plot_run = False

if plot_run == True:
    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)

    par1 = host.twinx()
    par2 = host.twinx()

    offset = 120
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right",
                    axes=par2,
                    offset=(offset, 0))

    par2.axis["right"].toggle(all=True)

    host.set_xlabel("Time [s]")
    host.set_ylabel(r"$K_{zx}$ [-]")
    par1.set_ylabel(r"$w_{zx}$ [1/s]")
    par2.set_ylabel(r"$z$ [-]")

    p1, = host.plot(data_arr[t_start_index:t_stop_index,0],data_arr[t_start_index:t_stop_index,4], 'r-', label=r'$u_{z}$')
    p2, = par1.plot(data_arr[t_start_index:t_stop_index,0],data_arr[t_start_index:t_stop_index,2], 'g-', label=r'$w_{zx}$')
    p3, = par2.plot(data_arr[t_start_index:t_stop_index,0],data_arr[t_start_index:t_stop_index,8], 'b-', label=r"z")

    host.legend(loc=3)

    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())
    par2.axis["right"].label.set_color(p3.get_color())

    matplotlib.rcParams.update({'font.size': 20})  # increase font size on axes (edited)
    plt.savefig('output/w_u_z_' + data_name + '.eps', bbox_inches='tight')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(data_arr[t_start_index:t_stop_index,0],data_arr[t_start_index:t_stop_index,9], 'b-', label=r"$x$")
    ax.set_ylabel("x [m]")
    ax.set_xlabel('Time [s]')
    plt.legend()
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    plt.savefig('output/x_' + data_name + '.eps', bbox_inches='tight')

    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.plot(data_arr[t_start_index:t_stop_index,0],data_arr[t_start_index:t_stop_index,5], 'b-', label=r"e cov($u_z,w_{zx}$)")
    ax.set_ylabel("covariance error [-]")
    ax.set_xlabel('Time [s]')
    plt.legend(loc=8)
    ax4 = ax.twinx()
    ax4.plot(data_arr[t_start_index:t_stop_index,0],data_arr[t_start_index:t_stop_index,7], 'g-', label=r'$K_{zx}$')
    ax4.set_ylabel(r'$K_{zx}$ [-]')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.legend(loc=4)
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    matplotlib.rcParams.update({'font.size': 20})  # increase font size on axes (edited)
    plt.savefig('output/cov_K_'+data_name+'.eps', bbox_inches='tight')

    plt.show()

else:
    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    ax.plot(day_2_x,day_2_K_avg, 'ro')
    ax.errorbar(day_2_x, day_2_K_avg, day_2_K_std, ls='none', color='red') #, color='black', elinewidth=2)
    ax.plot(x_range,K_range, 'g-')
    ax.set_ylabel(r'Instability gain $K_{zx}$ [-]')
    ax.set_xlabel('x distance [m]')
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    matplotlib.rcParams.update({'font.size': 20})  # increase font size on axes (edited)
    plt.savefig('output/linear_fake_flow.eps', bbox_inches='tight')

    plt.show()

