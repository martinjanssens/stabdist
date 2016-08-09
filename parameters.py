###############
#  PARAMETERS #
###############
import numpy as np

#MODES
#There are different modes. Setting the mode chooses the simulation
#Mode 1: Landing with constant gain, with delay, ZOH, constant wind, vz/z-based
#Mode 2: Hover until instability vz/z-based, with delay, constant wind, ZOH
#Mode 3: Hover until instability vz/x-based, with delay, constant wind, ZOH.
#        How are you going to deal with a camera looking horizontally when you
#        are moving up and down? How to assess when you become unstable,
#        especially when you are far away from the surface? Find out about OF-
#        methods in passing flow instead of divergence. Check linearity.
#Mode 4: Hover until instability vy/x-based, with delay, wind in y-direction, ZOH.
#        Better, worse or similar to vz/x? Still need to add the dynamics of tilting
#        the rotorcraft to generate uy!
#TO DO:
#   - Implement the dynamics of horizontal thrust and instability
#   - Rewrite "adaptgain" function so that gain is quickly reduced on detecting instability
#   - Add a symmetric error range for the covariance!

#Define mode
mode = 3

#Plot of a single run on or off. Set True if you run the simulate.py program
plot_traj = False
#Define whether you want to analyse (True) or not (False). Set True if running analysis.py
analysis = True
#Define whether you want to run distance.py (True) or not (False)
distance = False
#distance and analysis cannot both be True due to wind import in control function!

#General parameters
beta = 0.5    # --> beta = 1/2*rho*C_D*A
ti = 0.       #s
tlst = [ti]   #Initialise timelist
tmax = 55.    #s
g = 9.81      #m/s^2
PIinner = [1.,1.,-1.] #P,I,D gain for inner control loop 1.8,1.5
Iz = 1.         #Gain for integral control, inner control loop
covset = 0.01 #-
ecov = 0.005
emax = 0.0005 #Covariance error at which the program detects instability

#Variable initial parameters
delay = 0.1#s
Tstep = 0.01  #s
Kx0 = 30.      #- Initial gain
Kz0 = 10.      #- "
vwind = 0.    #m/s for single simulation
Kstep = 0.1  #s
tprev = ti    #s Previous time that the gain was adapted in outer loop

#Analysis of multiple hover manoeuvres with different initial conditions if analysis == True
if analysis == True or distance == True:
    #Analyse with different wind at different initial distances
    if mode == 2 or mode == 3 or mode == 4:  #only implemented in mode 2, 3 and 4
        # vwindlst = np.arange(-3.0,3.0,0.11)
        # x0lst = np.arange(8,12,2)
        # vwindlst = np.arange(-1.5,1.7,0.2)
        vwindlst = [0.]
        x0lst = np.arange(3,10,1)
        #Simulation function is based on a single value for wind and distance!
        
#Add gusts to the wind
gust = True
if gust == True:
    trange = np.arange(0,tmax+Tstep,Tstep)
    W = 0.1
    a = 3.0
    vwindgust = vwind + W*np.sin(a*trange)

#Noise model for divergence
noise = [0.8506,-0.0728,0.6022,0.1709,0.0306]

#STATE
#state = [x, vx, z, vz, ux, uz, m]
x = 10.
vx = 0.
z = 10.
vz = 0.
ux = 0.
uz = 0.
m = 1.        #kg
state = [x, vx, z, vz, ux, uz, m]

#State indices
ix = 0
ivx = 1
iz = 2
ivz = 3
iux = 4
iuz = 5
im = 6

#Initial optic flow (for perfect measurements)
wactx = vz/x
wactz = vz/z

#Out of range parameters
zmax = z+5.      #m
ymax = 50.      #m
ymin = 0.       #m
zland = 0.    #m

#Fourier analysis parameters
f_window = 35 #timesteps
bw_min = 1.5  #Hz
bw_max = 3.5  #Hz
Xset = 0.01#0.075
adaptcov = False
unstable = False

#Mode dependent parameters
if mode == 1:
    wset = -0.2     #1/s Set optic flow
           #[wactz, az, cov,ezlst]
    data = [[wactz],[g],[0.],[0.]]
    iwactz = 0
    iaz = 1
    icov = 2
    iezlst = 3
    PIouter = [0,0,0] #No outer loop
elif mode == 2:
    wset = 0.         #1/s
    PIouter = [.9,.2,.2] #P,I,D
          #[wactz,ecov,Kz,ez,tprev,Xfmaxlst,[wreal (without noise)]]
    data = [[wactz],[ecov],[Kz0],[0.],[tprev],[0.],[wactz]]
    iwactz = 0
    iecov = 1
    iKz = 2
    iezlst = 3
    itprev = 4
    iXfmax = 5
    iwrealz = 6
elif mode == 3:
    wset = 0.         #1/s
    PIouter = [.9,.2,.2] #P,I,D
          #[wactz,ecov,Kx,ez,tprev,drag,Xfmaxlst]
    data = [[wactx],[ecov],[Kx0],[0.],[tprev],[0.],[0.]]
    iwactx = 0
    iecov = 1
    iKx = 2
    iezlst = 3
    itprev = 4
    ip = 5
    iXfmax = 6
elif mode == 4:
    #Add the y-drection
    y = 10.
    vy = 0.
    iy = 7
    ivy = 8
    uy = 0.
    iuy = 9
    state.append(y)
    state.append(vy)
    state.append(uy)
    wacty = vy/x
    Ky = 10.
    wset = 0.         #1/s
    PIouter = [1,0,0] #P,I,D
          #[wacty,ecov,Ky,ey,tprev]
    data = [[wacty],[ecov],[Ky],[0.],[tprev]]
    iwacty = 0
    iecov = 1
    iKy = 2
    iezlst = 3 
    itprev = 4

#First state append
state_in_time = np.array(state)
