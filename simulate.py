################
#  SIMULATION  #
################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
import scipy.signal as scis
from scipy.integrate import odeint
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

def simulate(ti,tlst,Tstep,mode,state,data,state_in_time,tmax,zmax,zland,covset,vwind,delay,beta,Kstep,PIinner,PIouter,emax,f_window,Xset,adaptcov,noise):
    import control as c
    inrange = True
    while inrange:
        #Update time and define timespan for integration
        ti += Tstep
        tspan = [tlst[-1],ti]
        if mode == 1:
            from parameters import iuz,iz,icov
            #Control input and wind
            uz,p = c.controlmode1(state,data,ti,delay,Tstep,PIinner,beta,vwind)
            #Integration
            state = odeint(c.dydt_uzp,state,tspan,args=(uz,p,Tstep))[-1]
            state[iuz] = uz
            #Append the state and relevant data
            data, state_in_time = c.appendmode1(data,state,state_in_time,uz,p,noise)
            #Check for instability, landing or time runout
            if ti > tmax or state[iz] >= zmax or state[iz] <= zland or data[icov][-1] >= covset:
                inrange = False
        if mode == 2:
            from parameters import iuz,iz,iecov,iKz,iXfmax,unstable
            #Control law, wind, and K update
            uz,p,data = c.controlmode2(state,state_in_time,data,ti,vwind,delay,Tstep,beta,Kstep,covset,PIinner,PIouter,f_window,Xset,noise)
            #Integration
            state = odeint(c.dydt_uzp,state,tspan,args=(uz,p,Tstep))[-1]
            state[iuz] = uz
            #Append state and relevant data
            data, state_in_time = c.appendmode2(data,state,state_in_time,noise)
            #Check for instability, z out of range or time runout
            if adaptcov == False and data[iXfmax][-1]<=Xset-0.01 and len(data[iKz])>15:
                unstable = True
            elif adaptcov == True and (data[iecov][-1]<=covset-emax or data[iecov][-1]>=covset+emax) and len(data[iKz])>20:
                unstable = True
            if ti>tmax or state[iz]>=zmax or state[iz]<=zland or unstable == True:
                inrange = False
        if mode == 3:
            from parameters import iuz,iz,iecov,iKx,iXfmax,unstable
            #Control law, wind, Kz update
            uz,p,data = c.controlmode3(state,state_in_time,data,ti,vwind,delay,Tstep,beta,Kstep,covset,PIinner,PIouter,f_window,Xset,noise)
            #Integration
            state = odeint(c.dydt_uzp,state,tspan,args=(uz,p,Tstep))[-1]
            state[iuz] = uz
            #Append state and relevant data
            data, state_in_time = c.appendmode3(data,state,state_in_time,p,noise)
            #Check for instability after 10 K-updates (so initial wind isn't read as
            #instability), z out of range or time runout
            if adaptcov == False and data[iXfmax][-1]<=Xset-0.0005 and len(data[iKx])>25:
                unstable = True
            elif adaptcov == True and (data[iecov][-1]<=covset-emax or data[iecov][-1]>=covset+emax) and len(data[iKx])>20:
                unstable = True
            if ti>tmax or state[iz]>=zmax or state[iz]<=zland or unstable == True:
                inrange = False
        if mode == 4:
            from parameters import iuy,iy,iecov,iKy
            #Control law, wind, Kz update
            uy,py,Ky = c.controlmode4(state,state_in_time,data,ti,vwind,delay,Tstep,beta,Kstep,covset,PIinner,PIouter)
            #Integration
            state = odeint(c.dydt_uypy,state,tspan,args=(uy,py))[-1]
            state[iuy] = uy
            #Append state and relevant data
            data, state_in_time = c.appendmode4(data,state,state_in_time,noise)
            #Check for instability after 10 K-updates, y out of range or time runout
            if ti>tmax or state[iy]>=ymax or state[iy]<=ymin or (data[iecov][-1]<=emax and len(data[iKy])>10):
                inrange = False
        tlst.append(ti)
    return(data,state_in_time,tlst)

from parameters import plot_traj

if plot_traj == True: #Do a single simulation and plot the result
    from parameters import * #These are the starting parameters.
    data, state_in_time,tlst = simulate(ti,tlst,Tstep,mode,state,data,state_in_time,tmax,zmax,zland,covset,vwind,delay,beta,Kstep,PIinner,PIouter,emax,f_window,Xset,adaptcov,noise)
    if mode == 1:
        plt.subplot(1,3,1)
        plt.plot(tlst,state_in_time[:,iz])
        plt.title("z-t")
        plt.subplot(1,3,2)
        plt.plot(tlst,data[icov])
        plt.title("cov(uz,wz)")
        plt.subplot(1,3,3)
        plt.plot(tlst,data[iwactz])
        plt.title("vz/z")
        plt.show()
    if mode == 2:
        plt.subplot(1,4,1)
        plt.plot(tlst,state_in_time[:,iz])
        plt.title("z-t")
        plt.subplot(1,4,2)
        if adaptcov == True:
            plt.plot(data[itprev],data[iecov])
            plt.title("Covariance error")
        else:
            plt.plot(data[itprev],data[iXfmax])
            plt.title("|WACTX(f)|")
        plt.subplot(1,4,3)
        plt.plot(tlst,data[iwactz])
        plt.plot(tlst,data[iwrealz])
        plt.title("vz/z")
        plt.subplot(1,4,4)
        plt.plot(data[itprev],data[iKz])
        plt.title("Kx-t")
        plt.show()
    if mode == 3:
        # plt.subplot(1,2,1)
        # plt.plot(tlst,state_in_time[:,iz])
        # plt.title("z-t")
        # plt.subplot(1,4,2)
        # if adaptcov == True:
        #     plt.plot(data[itprev],data[iecov])
        #     plt.title("Covariance error")
        # else:
        #     plt.plot(data[itprev],data[iXfmax])
        #     plt.title("|WACTX(f)|")
        # plt.subplot(1,4,3)
        # plt.plot(tlst,data[iwactx])
        # plt.title("vz/x")
        # plt.subplot(1,4,4)
        # plt.plot(data[itprev],data[iKx])
        # plt.title("Kx-t")
        # plt.show()

        window_s = 15
        window = int(round(window_s/Tstep))+5
        window_k = int(round(window_s/Kstep))

        host = host_subplot(111, axes_class=AA.Axes)
        plt.subplots_adjust(right=0.75)

        par1 = host.twinx()
        par2 = host.twinx()
        # par3 = host.twinx()

        offset = 120
        new_fixed_axis = par2.get_grid_helper().new_fixed_axis
        par2.axis["right"] = new_fixed_axis(loc="right",
                                            axes=par2,
                                            offset=(offset, 0))

        par2.axis["right"].toggle(all=True)

        # new_fixed_axis = par3.get_grid_helper().new_fixed_axis
        # par3.axis["left"] = new_fixed_axis(loc="left",
        #                                     axes=par3,
        #                                     offset=(-offset, 0))
        #
        # par3.axis["left"].toggle(all=True)

        # host.set_xlim(0, 2)
        # host.set_ylim(0, 2)

        host.set_xlabel("Time")
        host.set_ylabel(r"$u_z$ [N]")
        par1.set_ylabel(r"$w_{zx}$ [1/s]")
        par2.set_ylabel(r"$e_{cov}$ [-]")
       # par3.set_ylabel(r"$z(t)$ [-]")

        p1, = host.plot(tlst[-window:], state_in_time[:, iuz][-window:], 'r-', label=r'$u_z$')
        p2, = par1.plot(tlst[-window:], data[iwactx][-window:], 'g-', label = r'$w_{zx}$')
        p3, = par2.plot(data[itprev][-window_k:], data[iecov][-window_k:], 'b-', label=r"$e_{cov}$")
        # p4, = host.plot(tlst[-window:], state_in_time[:, iz][-window:], 'k-', label=r'$z$')

        # par1.set_ylim(0, 4)
        # par2.set_ylim(bottom=0.0299)

        host.legend(loc=3)

        host.axis["left"].label.set_color(p1.get_color())
        par1.axis["right"].label.set_color(p2.get_color())
        par2.axis["right"].label.set_color(p3.get_color())
        # par3.axis["right"].label.set_color(p4.get_color())

        matplotlib.rcParams.update({'font.size': 20})  # increase font size on axes (edited)

        plt.draw()
        # plt.savefig('output\simulate_gust.eps', bbox_inches = 'tight')
        # plt.show()



        # window_s = 1.5
        # window = int(round(window_s/Tstep))+5
        # window_k = int(round(window_s/Kstep))
        # fig = plt.figure()
        # ax1 = fig.add_subplot(111)
        # # ax1.plot(tlst[-window:], state_in_time[:,iz][-window:], 'b-', label=r'$z(t)$')
        # ax1.plot(tlst[-window:], state_in_time[:, iuz][-window:], 'r-', label=r'$u_z$')
        # # ax1.plot(data[itprev][-window_k:], data[iXfmax][-window_k:], 'r-', label=r'$e_{|W|}$')
        # matplotlib.rcParams.update({'font.size': 20})  # increase font size on axes (edited)
        # # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        # # ax1.set_ylabel('$u_z$ [N]')
        # ax1.set_ylabel(r"Error in $|W_{ZX}(f)|$ [-]")
        # # ax1.set_ylim(top=0.031)
        # ax1.set_xlabel('Time [s]')
        # plt.legend(loc=3)
        # ax2 = ax1.twinx()
        # ax2.plot(tlst[-window:], data[iwactx][-window:], 'g-', label = r'$w_{zx}$')
        # # ax2.plot(tlst[-window:], state_in_time[:, ivz][-window:], 'b-', label=r'$v_z$')
        # # ax2.plot(data[itprev][-window_k:], data[iecov][-window_k:], 'b-', label=r"$e_{cov}$")
        # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        # ax2.set_ylabel(r'$w_{zx}$ [1/s], $v_z$ [m/s]')
        # plt.legend(loc=8)
        #
        # ax10 = ax1.twinx()
        # # ax10.plot(tlst[-window:], data[iwactx][-window:], 'g-', label=r'$w_{zx}$')
        # # ax2.plot(tlst[-window:], state_in_time[:, ivz][-window:], 'b-', label=r'$v_z$')
        # # ax2.plot(data[itprev][-window_k:], data[iecov][-window_k:], 'b-', label=r"$e_{cov}$")
        # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        # ax10.plot(data[itprev][-window_k:], data[iecov][-window_k:], 'b-', label=r"$e_{cov}$")
        # ax10.set_ylabel("Covariance error [-]")
        # plt.legend(loc=1)
        # plt.savefig('output\simulate_zoom_f.eps', bbox_inches = 'tight')

        fig2 = plt.figure()
        ax3 = fig2.add_subplot(111)
        if adaptcov == True:
            ax3.plot(data[itprev], data[iecov], 'b-', label=r"$e_{cov}$")
            ax3.set_ylabel("Covariance error [-]")
        else:
            ax3.plot(data[itprev], data[iXfmax], 'b-', label=r'$e_{|X|}$')
            ax3.set_ylabel(r"Error in $|W_{ZX}(f)|$ [-]") # Figure out the unit here!
        ax3.set_xlabel('Time [s]')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.legend(loc=8)
        ax4 = ax3.twinx()
        ax4.plot(data[itprev], data[iKx], 'g-', label=r'$K_{zx}$')
        ax4.set_ylabel(r'$K_{zx}$ [-]')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.legend(loc=4)
        # plt.savefig('output\simulate_e.eps', bbox_inches = 'tight')
        #
        # # Redefine the plotting window
        # window_s = 15.
        # window = int(round(window_s/Tstep))
        # window_k = int(round(window_s/Kstep))
        #
        # fig3 = plt.figure()
        # ax5 = fig3.add_subplot(111)
        # ax5.plot(tlst[-window:], state_in_time[:,iz][-window:], 'b-', label=r'$z(t)$')
        # # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        # ax5.set_ylabel('$z(t)$ [m]')
        # ax5.set_xlabel('Time [s]')
        # plt.legend(loc=8)
        # ax6 = ax5.twinx()
        # ax6.plot(tlst[-window:], data[iwactx][-window:], 'g-', label=r'$w_{zx}$')
        # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        # ax6.set_ylabel(r'$w_{zx}$ [1/s]')
        # plt.legend(loc=4)
        # plt.savefig('output\simulate_z.eps', bbox_inches='tight')

        plt.show()

    if mode == 4:
        plt.subplot(1,4,1)
        plt.plot(tlst,state_in_time[:,iy])
        plt.title("y-t")
        plt.subplot(1,4,2)
        plt.plot(data[itprev],data[iecov])
        plt.title("Error in covariance")
        plt.subplot(1,4,3)
        plt.plot(tlst,data[iwacty])
        plt.title("vy/x")
        plt.subplot(1,4,4)
        plt.plot(data[itprev],data[iKy])
        plt.title("Ky-t")
        plt.show()
            

#Detection of the divergence height depends heavily on where the instability
#is detected --> This is a cause for errors, but is reduced by quick detection.
#Problem doesn't exist when you estimate different axis
