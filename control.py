#############################
# Control functions and ode #
#############################

import numpy as np
import random

#Inner loop for constant OF with delay
def constOFctrl(ti,delay,wact,wset,wactlst,K,Tstep,ezlst,PIinner,Iz):
    if ti < delay:
        wact = wset  #Up until the delay time, the MAV senses wact as wset --> error = 0
        u = K*(wset - wact)
    else:
        e = wset - wactlst[-int(delay/Tstep)] #Then as the w from tdelay ago
        ezlst.append(e)
        up = PIinner[0]*K*e #Proportional control
        ui = PIinner[1]*Iz*sum(ezlst[:-int(delay/Tstep)]) #Integral control
        ud = PIinner[2]*(ezlst[-1]-ezlst[-2])/Tstep #Differential control
        u = ui + up + ud
    return u

#Calculate the OF without delay, but including noise
def OFestimator(v,pos,noise):
    from parameters import mode
    w = v/pos #Perfectly measured OF
    if mode == 2 or mode == 1 or mode == 3:
        #w_est = (w-noise[1])/noise[0]
        #w_est = w_est + noise[2]*w_est**2 + noise[3]*w_est# + noise[4]
        w_est = w + noise[2]*w**2 + noise[3]*w# Systematic noise
        w_est = w_est * float(random.randrange(98,103))/100
    elif mode == 4:
        #Random white noise, as divergence could be modelled like that too
        w_est = w*float(random.randrange(90,111))/100
    w_est = w # TURN OFF NOISE
    return w_est,w

#Outer loop for hover control
def adaptgain(wactlst,ulst,covset,K,PIouter,ti,ecovlst,delay,Tstep):
    #If the lists are empty or not long enough to work with the delay, set cov = 0
    if (len(wactlst) or len(ulst)) < int(delay/Tstep) + 10:
        cov = 0.0
    #Compute error in covariance at delay seconds ago
    else:
        ulst = np.array(ulst[-10:-1])
        wactlst = np.array(wactlst[-int(delay/Tstep)-10:-(int(delay/Tstep)+1)])
        #print ulst, wactlst
        #cov = np.cov(ulst,wactlst)[0,0] #cov is slightly positive for my delay setting!

        # Compute covariance
        uz_bar = np.average(ulst)
        wzx_bar = np.average(wactlst)

        diff_uz = ulst - uz_bar
        diff_wzx = wactlst - wzx_bar

        cov = sum(diff_uz * diff_wzx) / (ulst.shape[0] - 1)
        # print ti, cov
        # print "u: ",ulst
        # print "w: ",wactlst

    ecov = covset - cov   #A positive value in my case, it diverges negatively!
    # print "ecov: ", ecov
#Controllers to adapt the gain:
#I am struggling with building a controller that increases K enough to diverge, but not
#so much that it cannot re-adjust around the edge of instability
  #Option 1: What I understand from Guido's paper:
    K = K + PIouter[1]*K*ecovlst[-1] #K_current + Outergain_I*K_current*ecov_previous
    Kpi = K + PIouter[0]*K*ecov #K_new + Outergain_P*K_new*ecov_current
    Kpid = Kpi + PIouter[2]*Kpi*(ecov-ecovlst[-1]) #D-control
  #Option 2: PID controller --> Add differential gain to quicker adjust the gain
    #Kp = PIouter[0]*ecov
    #Ki = PIouter[1]*(sum(ecovlst)+ecov)
    #Kd = PIouter[2]*(ecov-ecovlst[-1])/Tstep #If negative, gain should be decreased
    #Kpi = K + Kp + Ki + Kd
    #Running does not show any significant improvement in ecov divergence after instability
  #Option 3: Just add a constant amount of gain!
    #Kpi = K + 1.9*(1+ecov)
    #Surprisingly gives the best results for a linear model, but will not work to
    #readjust the gain aroud the edge of instability
    if Kpid < 0:
        Kpid = 10.
    return Kpid, ecov

def adapt_K_freq(wactlst,K,PIouter,Tstep,delay,f_window,Xset,Xflst):
    from parameters import bw_min,bw_max
    #Fourier analysis
    if len(wactlst)-int(delay/Tstep) > f_window:
        x_window = wactlst[-int(delay/Tstep)-f_window:-int(delay/Tstep)]
        X_fft = abs(np.fft.rfft(x_window))
        f_X_fft = np.fft.rfftfreq(len(x_window), d=Tstep)
        #Pass everything in this range, nothing else
        X_fft[(f_X_fft<bw_min)] = 0.
        X_fft[(f_X_fft>bw_max)] = 0.
        #Track the max value in the designated range over time
        maxX = float(sum(X_fft))/len(X_fft)
    else:
        maxX = 0.
    #Use this signal for control:
    eXf = Xset - maxX #If positive, max is below threshold --> Icrease K
    K = K*(1 + PIouter[1]*Xflst[-1]) #K_current + Outergain_I*K_current*ecov_previous
    Kpi = K*(1 + PIouter[0]*eXf) #K_new + Outergain_P*K_new*ecov_current
    Kpid = Kpi*(1 + PIouter[2]*(eXf-Xflst[-1])) #D-control
    if Kpid < 0:
        Kpid = 10.
    return Kpid, eXf

def appendmode1(data,state,state_in_time,uz,p,noise):
    from parameters import iz,ivz,g,iwactz,iaz,icov,iuz,im,Kz
    wactz = OFestimator(state[ivz],state[iz],noise)
    az = (uz + p)/state[im] - g
    #If the lists are empty set cov = 0
    if len(data[iwactz]) < 2 or state_in_time.ndim < 2 :
        cov = 0.0
    else:
        cov = np.cov(data[iwactz],state_in_time[:,iuz])[0,0]*Kz/100.
    data[iwactz].append(wactz)
    data[iaz].append(az)
    data[icov].append(cov)
    state_in_time = np.vstack([state_in_time,state])
    return data, state_in_time

def controlmode1(state,data,ti,delay,Tstep,PIinner,beta,vwind):
    from parameters import gust,wset,Iz,g,Kz,iwactz,iezlst,im,ivz,iz,noise
    #Dynamics and control of z-axis
    wactz,_ = OFestimator(state[ivz],state[iz],noise)
    #Control law in z, including delay
    uz = state[im]*(g+constOFctrl(ti,delay,wactz,wset,data[iwactz],Kz,Tstep,data[iezlst],PIinner,Iz))
    #Drag
    if gust == False:
        p = np.sign(vwind-state[ivz])*0.5*beta*(vwind-state[ivz])**2
    if gust == True:
        from parameters import vwindgust
        p = np.sign(vwindgust[int(ti/Tstep)]-state[ivz])*0.5*beta*(vwindgust[int(ti/Tstep)]-state[ivz])**2
    return uz,p
    
def dydt_uzp (state,tspan,uz,p,Tstep):
    from parameters import g
    #Returns the derivave of the state vector and the control input,
    #for the case that uz is the input control variable

    #State indices, state = [x, vx, z, vz, ux, uz, m]
    ix = 0
    ivx = 1
    iz = 2
    ivz = 3
    iux = 4
    iuz = 5
    im = 6

    #Derivative of the state vector
    dstate_dt = np.zeros(7)
    dstate_dt[ix] = state[ivx]              #dx/dt = vx
    dstate_dt[ivx] = state[iux]/state[im]   #dvx/dt = ux/m <- no other force in x    
    dstate_dt[iuz] = uz - state[iuz]
    dstate_dt[iz] = state[ivz]
    dstate_dt[ivz] = (uz + p)/state[im] - g

    return dstate_dt
    
def controlmode2(state,state_in_time,data,ti,vwind,delay,Tstep,beta,Kstep,covset,PIinner,PIouter,f_window,Xset,noise):
    #Inner and outer loop control in mode 2
    from parameters import analysis,distance,gust,wset,Iz,g,iwactz,iezlst,im,ivz,iz,iuz,iecov,iKz,itprev,adaptcov,iXfmax
    #Dynamics and control of z-axis
    wactz,_ = OFestimator(state[ivz],state[iz],noise)
    #Inner Control law in z, including delay
    uz = state[im]*(g+constOFctrl(ti,delay,wactz,wset,data[iwactz],data[iKz][-1],Tstep,data[iezlst],PIinner,Iz))
    #Outer control law in z if more that Kstep seconds have passed since previous update
    if ti > data[itprev][-1] + Kstep:
        if adaptcov == True:
            Kz,ecov = adaptgain(data[iwactz][1:],state_in_time[1:,iuz]-g,covset,data[iKz][-1],PIouter,ti,data[iecov],delay,Tstep)
            data[iecov].append(ecov)
        else:
            Kz,eXf = adapt_K_freq(data[iwactz],data[iKz][-1],PIouter,Tstep,delay,f_window,Xset,data[iXfmax])
            data[iXfmax].append(eXf)
        tprev = ti
        data[iKz].append(Kz)
        data[itprev].append(tprev)
    #Drag, including gusts. Import actual gust value if analysing different conditions
    if gust == False:
        p = np.sign(vwind-state[ivz])*0.5*beta*(vwind-state[ivz])**2
    elif gust == True and analysis == False:
        from parameters import vwindgust
        p = np.sign(vwindgust[int(ti/Tstep)]-state[ivz])*0.5*beta*(vwindgust[int(ti/Tstep)]-state[ivz])**2
    elif gust == True and analysis == True:
        from analysis import vwindgust
        p = np.sign(vwindgust[int(ti/Tstep)]-state[ivz])*0.5*beta*(vwindgust[int(ti/Tstep)]-state[ivz])**2
    elif gust == True and distance == True:
        from distance import vwindgust
        p = np.sign(vwindgust[int(ti/Tstep)]-state[ivz])*0.5*beta*(vwindgust[int(ti/Tstep)]-state[ivz])**2
    return uz,p,data

def appendmode2(data,state,state_in_time,noise):
    #Appends the state and wz in mode 2
    from parameters import iz,ivz,iwactz,iwrealz
    wactz,wreal = OFestimator(state[ivz],state[iz],noise)
    data[iwactz].append(wactz)
    data[iwrealz].append(wreal)
    state_in_time = np.vstack([state_in_time,state])
    return data, state_in_time

def controlmode3(state,state_in_time,data,ti,vwind,delay,Tstep,beta,Kstep,covset,PIinner,PIouter,f_window,Xset,noise):
    #Inner and outer loop control in mode 3
    from parameters import adaptcov,analysis,gust,iezlst,wset,Iz,g,im,ivz,ix,iuz,iecov,iKx,itprev,iwactx,distance,iXfmax
    #Dynamics and control of z-axis, given perfect knowledge of vz/x
    wactx,_ = OFestimator(state[ivz],state[ix],noise)
    #Inner Control law in z, including delay
    uz = constOFctrl(ti,delay,wactx,wset,data[iwactx],data[iKx][-1],Tstep,data[iezlst],PIinner,Iz)
    #Scale to thrust
    uz = state[im] * (g + uz)
    #Add actuator efficiency
    uz_eff = uz - 0.5*state[ivz]*uz - 0.5*state[ivz]
    if uz_eff > 0:
        uz = uz_eff

    #Add saturation
    u_sat = 10.4
    if uz > u_sat: uz = u_sat
    elif uz < -u_sat: uz = -u_sat

    #Outer control law in z
    if ti > data[itprev][-1] + Kstep:
        if adaptcov == True:
            Kx,ecov = adaptgain(data[iwactx][1:],state_in_time[1:,iuz]-g,covset,data[iKx][-1],PIouter,ti,data[iecov],delay,Tstep)
            data[iecov].append(ecov)
        else:
            Kx,eXf = adapt_K_freq(data[iwactx],data[iKx][-1],PIouter,Tstep,delay,f_window,Xset,data[iXfmax])
            data[iXfmax].append(eXf)
        tprev = ti
        data[iKx].append(Kx)
        data[itprev].append(tprev)
    #Drag, including gusts. Import actual gust value if analysing different conditions
    if gust == False:
        p = np.sign(vwind-state[ivz])*0.5*beta*(vwind-state[ivz])**2
    elif gust == True and analysis == False:
        from parameters import vwindgust
        p = np.sign(vwindgust[int(ti/Tstep)]-state[ivz])*0.5*beta*(vwindgust[int(ti/Tstep)]-state[ivz])**2
    elif gust == True and analysis == True:
        from analysis import vwindgust
        p = np.sign(vwindgust[int(ti/Tstep)]-state[ivz])*0.5*beta*(vwindgust[int(ti/Tstep)]-state[ivz])**2
    elif gust == True and distance == True:
        from distance import vwindgust
        p = np.sign(vwindgust[int(ti/Tstep)]-state[ivz])*0.5*beta*(vwindgust[int(ti/Tstep)]-state[ivz])**2
    return uz,p,data

def appendmode3(data,state,state_in_time,p,noise):
    #Appends the state and wx in mode 3
    from parameters import ix,ivz,iwactx,ip,f_window,bw_min,bw_max,Tstep,iXfmax
    wactx,_ = OFestimator(state[ivz],state[ix],noise)
    data[iwactx].append(wactx)
    data[ip].append(p)
    state_in_time = np.vstack([state_in_time,state])
    return data, state_in_time

def controlmode4(state,state_in_time,data,ti,vwind,delay,Tstep,beta,Kstep,covset,PIinner,PIouter):
    #Inner and outer loop control in mode 3
    from parameters import analysis,gust,iezlst,wset,Iz,im,ix,iuy,ivy,iecov,iKy,itprev,iwacty,distance
    #Dynamics and control of y-axis, including noise
    wacty,_ = OFestimator(state[ivy],state[ix],noise)
    #Inner Control law in y, including delay
    uy = state[im]*constOFctrl(ti,delay,wacty,wset,data[iwacty],data[iKy][-1],Tstep,data[iezlst],PIinner,Iz)
    #Outer control law in y
    if ti > data[itprev][-1] + Kstep:
        Ky,ecov = adaptgain(data[iwacty][1:],state_in_time[1:,iuy],covset,data[iKy][-1],PIouter,ti,data[iecov],delay,Tstep)
        tprev = ti
        data[iecov].append(ecov)
        data[iKy].append(Ky)
        data[itprev].append(tprev)
    #Drag, including gusts (same gusts as in z). Import actual gust value if analysing different conditions
    if gust == False:
        p = np.sign(vwind-state[ivy])*0.5*beta*(vwind-state[ivy])**2
    elif gust == True and analysis == False:
        from parameters import vwindgust
        p = np.sign(vwindgust[int(ti/Tstep)]-state[ivy])*0.5*beta*(vwindgust[int(ti/Tstep)]-state[ivy])**2
    elif gust == True and analysis == True:
        from analysis import vwindgust
        p = np.sign(vwindgust[int(ti/Tstep)]-state[ivy])*0.5*beta*(vwindgust[int(ti/Tstep)]-state[ivy])**2
    elif gust == True and distance == True:
        from distance import vwindgust
        p = np.sign(vwindgust[int(ti/Tstep)]-state[ivz])*0.5*beta*(vwindgust[int(ti/Tstep)]-state[ivz])**2
    return uy,p,data

def dydt_uypy (state,tspan,uy,py):
    from parameters import g
    #Returns the derivave of the state vector and the control input

    #State indices, state = [x, vx, z, vz, ux, uz, m, y, vy, uy]
    ix = 0
    ivx = 1
    iz = 2
    ivz = 3
    iux = 4
    iuz = 5
    im = 6
    iy = 7
    ivy = 8
    iuy = 9

    #Derivative of the state vector
    dstate_dt = np.zeros(10)
    dstate_dt[ix] = state[ivx]              #dx/dt = vx
    dstate_dt[ivx] = state[iux]/state[im]   #dvx/dt = ux/m <- no other force in x
    dstate_dt[iy] = state[ivy]              #dy/dt = vy
    dstate_dt[iuy] = uy - state[iuy]
    dstate_dt[ivy] = (uy + py)/state[im]    #No other acceleration in y
    uz = g
    dstate_dt[iuz] = uz - state[iuz]
    dstate_dt[iz] = state[ivz]
    dstate_dt[ivz] = state[iuz]/state[im] - g

    return dstate_dt

def appendmode4(data,state,state_in_time,noise):
    #Appends the state and wx in mode 3
    from parameters import ix,ivy,iwacty
    wacty,_ = OFestimator(state[ivy],state[ix],noise)
    data[iwacty].append(wacty)
    state_in_time = np.vstack([state_in_time,state])
    return data, state_in_time
