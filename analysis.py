##############
#  ANALYSIS  #
##############

import numpy as np
import matplotlib.pyplot as plt
import math
from simulate import simulate
import scipy.stats as sst
import sys
import matplotlib

from parameters import analysis,mode,x0lst,vwindlst

Kcritlst = []
xcritlst = []
if analysis == True:
    #Only implemented for mode 2, 3 and 4
    if mode == 2 or mode == 3 or mode == 4:
        for x0 in x0lst:
            for v in vwindlst:
                #Reset all the parameters to initial parameters
                exec(open('parameters.py').read())
                #Define wind and distance settings
                vwind = v
                vwindgust = vwind + W*np.sin(a*trange)
                if mode == 3 or mode == 4:
                    state[ix] = x0
                elif mode == 2:
                    state[iz] = x0
                #Do the simulation
                data,state_in_time,tlst = simulate(ti,tlst,Tstep,mode,state,data,state_in_time,tmax,zmax,zland,covset,vwind,delay,beta,Kstep,PIinner,PIouter,emax,f_window,Xset,adaptcov,noise)
                #Append Kx and x on instability point, at ti - delay
                if mode == 2:
                    Kcritlst.append(data[iKz][-int(delay/Tstep)])
                    xcritlst.append(state_in_time[:,iz][-int(delay/Tstep)])
                    print "Done: z0 = ",x0,"vwind = ",vwind,". Kz = ",data[iKz][-int(delay/Tstep)]
                elif mode == 3:
                    Kcritlst.append(data[iKx][-int(delay/Tstep)])
                    xcritlst.append(state_in_time[:,ix][-int(delay/Tstep)])
                    print "Done: x0 = ",x0,"vwind = ",vwind,". Kx = ",data[iKx][-int(delay/Tstep)] 
                elif mode == 4:
                    Kcritlst.append(data[iKy][-int(delay/Tstep)])
                    xcritlst.append(state_in_time[:,ix][-int(delay/Tstep)])
                    print "Done: x0 = ",x0,"vwind = ",vwind,". Kx = ",data[iKy][[-int(delay/Tstep)]]
                #Check for the effect of wind on the theoretically predicted slope
                #ppeak = max(abs(max(data[ip][:int(2/Tstep)])),abs(min(data[ip][:int(2/Tstep)])))
                #if ppeak != 0.:
                #    print(math.exp(ppeak*Tstep)-1)/(ppeak*(math.exp(ppeak*Tstep)+1)), 0.5*Tstep
    #Statistical analysis
    slope,intercept,r,pvalue,std = sst.linregress(Kcritlst,xcritlst)
    print "slope:", slope
    print "intercept:", intercept
    print "r^2:", r**2
    print "p-value:", pvalue
    if len(xcritlst) < len(x0lst)*len(vwindlst):
        print "Time runouts occurred!"

    #Plot a line of best fit
    x_bestfit = np.arange(min(Kcritlst),max(Kcritlst),0.1)
    y_bestfit = intercept + slope*x_bestfit

    #Theoretical distance prediction: x = 0.5*Kx*Tstep
    #ppeak = max(data[ip][:int(2/Tstep)]) #Pick out the first peak
    #Kxtheoretical = np.arange(0.,max(Kcritlst),0.1)
    #xcrittheoretical = Kxtheoretical*Tstep*0.5
    #coeff_p = (math.exp(ppeak*Tstep)-1)/(ppeak*(math.exp(ppeak*Tstep)+1))
    #xcrit_p = Kxtheoretical*coeff_p

    #Plot
    plt.plot(Kcritlst,xcritlst, "ro")
    # fig1 = plt.figure()
    # ax = fig1.add_subplot(111)
    # for x0 in range(len(x0lst)):
    #     xcrit = xcritlst[x0*len(vwindlst):x0*len(vwindlst)+len(vwindlst)]
    #     Kcrit = Kcritlst[x0*len(vwindlst):x0*len(vwindlst)+len(vwindlst)]
    #     plt.scatter(Kcrit, xcrit, c=vwindlst, marker='o', s=50, cmap='winter')
    plt.plot(x_bestfit,y_bestfit, "g-")
    plt.xlabel(r"$K_{zx}$ [-]")
    plt.ylabel(r"$x$ [m]")
    matplotlib.rcParams.update({'font.size': 20})  # increase font size on axes (edited)
    # plt.colorbar(label='Wind speed [m/s]')

    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    #plt.plot(Kxtheoretical,xcrittheoretical, "b-")
    #plt.plot(Kxtheoretical,xcrittheoretical, "m-")
    plt.savefig('output\dist_analysis_f.eps', bbox_inches='tight')
    plt.show()

    # Plot to show relation with wind
    # for x0 in range(len(x0lst)):
    #     Kcrit = Kcritlst[x0 * len(vwindlst):x0 * len(vwindlst) + len(vwindlst)]
    #     plt.scatter(vwindlst, Kcrit)
    #
    # plt.show()

    #Estimate the distance at a further point based on the compiled model, no wind
    x_real = x0lst[-1] + 2.
    exec(open('parameters.py').read())
    vwind = 0.5
    vwindgust = vwind + W*np.sin(a*trange)
    if mode == 3 or mode == 4:
        state[ix] = x_real
    if mode == 2:
        state[iz] = x_real
    data,state_in_time,tlst = simulate(ti,tlst,Tstep,mode,state,data,state_in_time,tmax,zmax,zland,covset,vwind,delay,beta,Kstep,PIinner,PIouter,emax,f_window,Xset,adaptcov,noise)
    if mode == 2:
        Kcrit = data[iKz][-int(delay/Tstep)]
    if mode == 3:
        Kcrit = data[iKx][-int(delay/Tstep)]
    elif mode == 4:
        Kcrit = data[iKy][-int(delay/Tstep)]
    x_predicted = intercept + slope*Kcrit
    print x_real, x_predicted
sys.exit()
