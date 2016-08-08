########################
#  DISTANCE PREDICTOR  #
########################

import numpy as np
import matplotlib.pyplot as plt
import math
from simulate import simulate
import scipy.stats as sst

#Define whether you want to run the parameter effect test or the model estimation
distance_pareffect = True
distance_modelest = False

#Test of the different parameters' effect on the model

#If I can show that the linear Kx-x model is purely a function of things that I can
#control and pre-program (and not of external variables like wind), then you have
#a very good distance estimator!
#What are the inputs? Accuracy of the OF measurement, covariance setting,
#mass of the drone, gain settings of the update loops, beta, delay, timestep,
#how often the gain is updated

from parameters import *
#Only run if so commanded from parameters.py
if distance == True and distance_pareffect == True:
    #Define the parameters I want to test
    testpars = [covset,m,Kstep,Tstep,delay,beta,PIinner[0],PIinner[1],PIouter[0],PIouter[1]]
    covset_i = 0
    m_i = 1
    Kstep_i = 2
    Tstep_i = 3
    delay_i = 4
    beta_i = 5
    PIinner0_i = 6
    PIinner1_i = 7
    PIouter0_i = 8
    PIouter1_i = 9
    
    outputs = np.zeros(5)
    for index in range(len(testpars)):       
        #Test each parameter with a  20% decrease and increase
        testlst = [testpars[index]*0.8,testpars[index],testpars[index]*1.2]
        for test in range(len(testlst)):
            Kcritlst = []
            xcritlst = []
            for x0 in x0lst:
                for v in vwindlst:
                    exec(open('parameters.py').read())                
                    #Define wind and distance settings
                    vwind = v
                    vwindgust = vwind + W*np.sin(a*trange)
                    state[ix] = x0
                    #Define the different tests (need ifs to overrule parameters.py import)
                    if index == covset_i:
                        covset = testlst[test]
                        ecov = covset
                        emax = covset - 0.001
                    elif index == m_i:
                        state[im] = testlst[test]
                    elif index == Kstep_i:
                        Kstep = testlst[test] 
                    elif index == Tstep_i:
                        Tstep = testlst[test]
                    elif index == delay_i:
                        delay = testlst[test] 
                    elif index == beta_i:
                        beta = testlst[test]
                    elif index == PIinner0_i:
                        PIinner[0] = testlst[test]
                    elif index == PIinner1_i:
                        PIinner[1] = testlst[test]
                    elif index == PIouter0_i:
                        PIouter[0] = testlst[test]
                    elif index == PIouter1_i:
                        PIouter[1] = testlst[test]
                    #Do the simulation
                    data,state_in_time,tlst = simulate(ti,tlst,Tstep,mode,state,data,state_in_time,tmax,zmax,zland,covset,vwind,delay,beta,Kstep,PIinner,PIouter,emax,f_window,noise)
                    #Append K and x on instability point
                    if mode == 3:
                        Kcritlst.append(data[iKx][-int(delay/Tstep)])
                    elif mode == 4:
                        Kcritlst.append(data[iKy][-int(delay/Tstep)])
                    xcritlst.append(state_in_time[-int(delay/Tstep)][ix])
            #Make model
            slope,intercept,r,pvalue,std = sst.linregress(Kcritlst,xcritlst)
            #Add the models and corresponding test
                    #[testparameter, value of testparameter, slope, intercept, std]
            output = [index, testlst[test], slope, intercept, std]
            outputs = np.vstack([outputs,output])
            print "Done: index = ",index,"test = ",testlst[test],"Percent:",100.*(len(testlst)*index+test+1)/(len(testpars)*len(testlst)),"%"

#Results of running with mode = 3   [covset,m, Kstep,Tstep,delay,beta]
#Initial parameters:                [0.03, 1.0, 0.15, 0.03, 0.15, 0.5]
#        index, value, slope,           intercept,      r^2      
##np.array([[0.0, 0.024, 0.0472420587098, 0.375251038248, 0.999199296323], \\
##          [0.0, 0.03,  0.0451426809392, 0.319158577526, 0.999365648103], \\
##          [0.0, 0.036, 0.042384321904,  0.182357797173, 0.998353219961], \\
##          [1.0, 0.8,   0.0434160728032, 0.30153199767,  0.99827205666], \\
##          [1.0, 1.0,   0.0434302572365, 0.301584722466, 0.998277091661], \\
##          [1.0, 1.2,   0.0434398597018, 0.301630855656, 0.99828068976], \\
##          [2.0, 0.12,  0.0434476033354, 0.054149695902, 0.999751310487], \\
##          [2.0, 0.15,  0.0434302572365, 0.301584722466, 0.998277091661], \\
##          [2.0, 0.18,  0.0442545957652, 0.387322669256, 0.998433015076], \\
##          [3.0, 0.024, 0.0439364192331, 0.247678356651, 0.999708872024], \\
##          [3.0, 0.03,  0.0434302572365, 0.301584722466, 0.998277091661], \\
##          [3.0, 0.036, 0.0415682706607, 0.272896968897, 0.999610715679], \\
##          [4.0, 0.12,  0.0354547184596, 0.297104526231, 0.997679910458], \\
##          [4.0, 0.15,  0.0434302572365, 0.301584722466, 0.998277091661], \\
##          [4.0, 0.18,  0.0509315329564, 0.380009476322, 0.999127365999], \\
##          [5.0, 0.4,   0.0430791513048, 0.338732831067, 0.998355742503], \\
##          [5.0, 0.5,   0.0434302572365, 0.301584722466, 0.998277091661], \\
##          [5.0, 0.6,   0.0436544560759, 0.300474352034, 0.998136539157]])
#Indicates that :
#   - Increasing covset: Decreases the slope, decreases the y-intercept
#   - Increasing m:      NO EFFECT
#   - Increasing Kstep:  NO EFFECT on slope; change in y-intercept; worse accuracy
#   - Increasing Tstep:  Slightly decreases the slope, decreases y-intercept
#   - Increasing delay:  Increases slope, increases y-intercept
#   - Increasing beta:   NO EFFECT
#   - Increasing inner P-gain: Increases slope, decreases y-intercept, increases accuracy
#   - Increasing inner I-gain: No effect on slope, increases y-intercept
#   - Increasing outer P-gain: No effect on slope, decreases y-intercept

#Which means that for a certain drone with a given setting of [covset,Kstep,Tstep,delay],
#you should be able to implement a linear distance estimator that works independently of
#the wind-conditions, the mass of the drone (so with payload!), air density, surface area
#or drag! That's nice!


#Perform two simulations at known distances, and make a model.
#Based on the model perform another simulation and use it to estimate the distance
#based purely on that model. Not very useful right now...
elif distance == True and distance_modelest == True:
    #With random wind settings perform two simulations at known distances, and
    #make a model. Based on the model perform another simulation and use it to
    #estimate the distance based purely on the model.
    Kcritlst = []
    xcritlst = []
    #How to choose at which heights to make the estimates?
    xtest = [2.,6.]
    for x0 in xtest:
        exec(open('parameters.py').read())
        #Define wind and distance settings
        vwind = 0.
        #vwindgust = vwind + W*np.sin(a*trange)
        state[ix] = x0
        #Do the simulation
        data,state_in_time,tlst = simulate(ti,tlst,Tstep,mode,state,data,state_in_time,tmax,zmax,zland,covset,vwind,delay,beta,Kstep,PIinner,PIouter,emax,f_window,Xset,adaptcov,noise)
        #Append Kx and x on instability point
        if mode == 3:
            Kcritlst.append(data[iKx][-int(delay/Tstep)])
            xcritlst.append(state_in_time[:,ix][-int(delay/Tstep)])
            print "Done: x0 = ",x0,"vwind = ",vwind,". Kx = ",data[iKx][-int(delay/Tstep)]
        elif mode == 4:
            Kcritlst.append(data[iKy][-int(delay/Tstep)])
            xcritlst.append(state_in_time[:,ix][-int(delay/Tstep)])
            print "Done: x0 = ",x0,"vwind = ",vwind,". Kx = ",data[iKy][[-int(delay/Tstep)]]
    #Model based on the runs
    slope,intercept,r,pvalue,std = sst.linregress(Kcritlst,xcritlst)
    #Try the model
    x_real = 10.
    exec(open('parameters.py').read())
    vwind = 0.
    vwindgust = 0.
    state[ix] = x_real
    data,state_in_time,tlst = simulate(ti,tlst,Tstep,mode,state,data,state_in_time,tmax,zmax,zland,covset,vwind,delay,beta,Kstep,PIinner,PIouter,emax,f_window,Xset,adaptcov)
    if mode == 3:
        Kcrit = data[iKx][-int(delay/Tstep)]
    elif mode == 4:
        Kcrit = data[iKy][-int(delay/Tstep)]
    x_predicted = intercept + slope*Kcrit
    print x_real, x_predicted
    
#What I may want to know:
#  - The effect of increasing the distance between the two points
#  - The effect of increasing the number of test runs before compiling a model

