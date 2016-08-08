# stabdist
Simulation environment for a stability-based, monocular distance estimation algorithm for UAV systems 

MODES <br />
There are different modes, set in parameters.py. Setting the mode chooses the simulation. <br />
Mode 1: Landing with constant gain, with delay, ZOH, constant wind, vz/z-based <br />
Mode 2: Hover until instability vz/z-based, with delay, constant wind, ZOH. <br />
Mode 3: Hover until instability vz/x-based, with delay, constant wind, ZOH. <br />
Mode 4: Hover until instability vy/x-based, with delay, wind in y-direction, ZOH. <br />

Set plot_traj = True and analysis = False to run/plot a single simulation run. <br />
Conversely set analysis = True and plot_traj = False to run a set of simulations for a range of distances and wind speeds.
