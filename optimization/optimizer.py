# could improve:
# minimizer
# actual analysis

# Nelder-Mead simplex search

import numpy as np
import trajectory
import openrocket_interface as openrkt
from math import sqrt, pi, exp, log, cos
import math as m # is this sanitary? is this needed? 'm' is a variable name too...
from scipy.optimize import minimize
global cons_TWR, cons_S_crit, cons_accel, cons_LD, cons_alt, X0, m0 # i'm sorry for global vars...
global allvectors, cons_mass, dbz
dbz = 0

# SIMULATION AND OPTIMIZATION PARAMETERS
time_step = 0.1 # change time-step for trajectory simulation
iterations = 10 # number of escalating iterations, degenerates after ~8
launch_site_alt = 1401 # m, altitude of launch site above sea level

##CHANGE INITIAL DESIGN GUESS HERE
# be sure that you start with a feasible design, otherwise the problem will be ill-conditioned
L = 1.7    # Total tank lengths (m)
mdot = 2.2 # Propellant mass flow rate (kg/s)
dia = 12.  # Rocket diameter (in)
p_e = 47.  # Exit Pressure (kPa)

#CHANGE CONSTRAINTS HERE
cons_mass = 200.                    # GLOW constraint, kg
cons_ls = 22.                       # min launch speed from 60' tower constraint, m/s
cons_TWR = 2.                        # TWR constraint
cons_S_crit = 0.35                  # Critical pressure ratio constraint
cons_accel = 15.                    # Max acceleration constraint, g's
cons_LD = 18.                       # L/D ratio constraint
cons_alt = 100000. + launch_site_alt # Min altitude constraint, km
cons_thrust = 6.                    # max ground-level thrust, kN
base11_ceiling = 150000.             # base-11 maximum apogee requirement, km

shirt_l = trajectory.shirt_l # non-engine subsystems lengths
g_0 = openrkt.g_0 # gravity

X0 = np.array([L, mdot, p_e]) #numpy arrays are nicer
allvectors = [] # array for all design vecs

# calculates length-diameter ratio
def ld_ratio(pants_l, dia):
    total_l = shirt_l + pants_l # adds shirt and pants lengths to get total rocket length
    D = dia * 0.0254 # converts in to m for total airframe diameter
    return total_l / D # gets ratio
    
# calculates top g-force
def max_g_force(a):
    max_gees = max(abs(a)) / g_0
    return max_gees
    
#gets the index of the last thrust term before engine kicks
def last_moment(F):
    j = 0
    fdex = 1 # init fdex is a kludge because optimizer crashed once
    
    for thing in F:
        if thing == 0:
            fdex = j
            break
        j += 1
    return fdex

# all of our comparisons are ratios instead of subtractions because
# it's normalized, instead of sensitive to magnitudes of constraints

# minimize this, **2 makes it convex and well behaved w.r.t. when var=cons
def objective(var, cons):
    return (var/cons)**2 / 2

# **2 for convexity and because i like it more than abs(), but that also works
def exact(var, cons):
    return (var/cons - 1)**2 / 2

# this is your basic exterior penalty, either punishes for unfeasibility or is inactive
def exterior(var, cons, sat_if_less_than=False):
    if sat_if_less_than:
        return max(0, var/cons - 1)**2 / 2
    else:
        return max(0, -(var/cons - 1))**2 / 2

# this restricts our objective function to the strictly feasible region
# make rockets great again, build that wall, etc
# i like the logorithmic version more, but interior point method is viable too
def barrier(var, cons, int_point=False, sat_if_less_than=True):
    global dbz
    def interior(g):
        return 1/g
    try:
        #print('log: '+str(-log((var/cons - 1)))+ ', 1/x: '+str(-interior(-(var/cons - 1)))) # print for comparison
        if not int_point:
            if sat_if_less_than:
                #if var > cons:
                    #print "division by zero!"
                return -log(-(var/cons - 1))
            else:
                #if var < cons:
                    #print "division by zero!"
                return -log((var/cons - 1))
        else:
            if sat_if_less_than:
                return -interior(var/cons - 1)
            else:
                return -interior(-(var/cons - 1))
    except:
        dbz += 1
        return 10**8

# this manages all our constraints, penalty parameters: mu zooms out and rho zooms in 
def penalty(ls, F, LD, TWR, S_crit, alt, max_g, mu, rho):
    b = [barrier(alt, cons_alt,False,False), barrier(alt, base11_ceiling,False, True) ]
    eq = [exact(F, cons_thrust)]
    ext = [exterior(ls, cons_ls) , exterior(LD, cons_LD, True) , exterior(TWR, cons_TWR) , \
          exterior(S_crit, cons_S_crit) , exterior(max_g, cons_accel, True)]# , \
          #exterior(alt, cons_alt)]
    #print("barrier: " + str(mu*sum(b))+", exact: "+ str(rho/2 *sum(eq))+", exterior: "+ str(rho*sum(ext)))
    return mu*sum(b) + rho*(sum(eq)/2 + sum(ext))

# Pseudo-objective function
# x is array of design parameters, n is sequence index of penalty and barrier functions
# print blocks are sanity checks so i'm not staring at a blank screen and can see what various tweaks actually do
def f(x, n):
    global allvectors
    L = x[0]   # Tank length (m)
    mdot = x[1] # Propellant mass flow rate (kg/s)
    p_e = x[2] #x[3]  # Pressure (kPa)
    
    # get trajectory data from x
    sim = trajectory.trajectory(L, mdot, dia, p_e, x_init=launch_site_alt, dt=time_step)
    fdex = last_moment(sim.F) # index of engine burnout
    r, l_o, l_f = openrkt.split_tanks(sim.m_prop[0], dia)
    obj = objective(sim.m[0], cons_mass) # minimize GLOW
    #print("obj: "+str(obj))
    # then, calculate penalization from trajectory based on initial thrust
    pen = penalty(sim.launch_speed, sim.F[0]/1000, ld_ratio(openrkt.system_length(l_o, l_f), dia), \
        sim.TWR, sim.S_crit, sim.alt[-1], max_g_force(sim.a),(.8)/(2**n),(1.25)*(2**n)) # mu and rho are selected for nice behavior
    
    # add objective and penalty functions
    sum_func = obj + pen
    #print("Total (n): "+str(sum_func) + ' ('+ str(n) +')')
    #print("L "+ str(L)+ ", D " + str(dia)+ ", mdot " + str(mdot)+ ", p_e "+  str(p_e)+ ", alt "+  str(sim.alt[-1]))
    #print("obj "+ str(obj)+ ' + pen '+ str(pen)+ ' = '+ str(sum_func))
    #print('')
    allvectors.append(x) # maintains a list of every design
    return sum_func

# we want to iterate our optimizer for convergence reasons  
# n = number of sequential iterations 
def iterate(f, x_0, n):
    x = x_0 # initial design vector
    global dbz
    for i in range(n):
        print("Iteration " + str(i+1) + ":")
        res = minimize(f, x, args=(i+0), method='nelder-mead', options={'disp': True, 'adaptive':True}) # this minimizer uses simplex method
        x = res.x # feed optimal design vec into next iteration
        alt = trajectory.trajectory(x[0], x[1], dia, x[2], x_init=launch_site_alt, dt=time_step).alt
        print("         Divisions by zero (extreme violations of altitude window): "+str(dbz))
        print("Propellant tube length (m): "+str(x[0]))
        print("Mass flow rate (kg/s): "+str(x[1]))
        print("Airframe diameter (in): "+str(12))
        print("Exit pressure (kPa): "+str(x[2]))
        print("Altitude (km): "+str(alt[-1]))
        print('')
        dbz=0
    return x

# this creates a list of strings for relevant data of trajectory
def print_results(res):
    text_base = [] # list of lines of strings
    
    # Rename the optimized output
    L = res[0]
    mdot = res[1]
    p_e = res[2]
    
    sim = trajectory.trajectory(L, mdot, dia, p_e, x_init=launch_site_alt, dt=time_step) 
    fdex = last_moment(sim.F) # index of engine burnout
    r, l_o, l_f = openrkt.split_tanks(sim.m_prop[0], dia)
    eng_sys_len = openrkt.system_length(l_o, l_f)
    
    np.set_printoptions(precision=3) # this line may be deprecated, i copy-pasted most of this section
    
    text_base.append('\nOPTIMIZED DESIGN VECTOR')
    text_base.append('\n-----------------------------')
    text_base.append('\nx_initial_guess                            = ' + ', '.join([str(X0[0]), str(X0[1]), str(dia), str(X0[2])])) # kind of a kludge
    text_base.append('\ninitial guess GLOW                                    = {:.1f} kg'.format( \
          trajectory.trajectory(X0[0], X0[1], dia, X0[2], x_init=launch_site_alt, dt=time_step).m[0]))
    text_base.append('\nx_optimized                                = ' + ', '.join([str(L), str(mdot), str(dia), str(p_e)])) # kind of a kludge
    text_base.append('\ndesign GLOW                                = {:.1f} kg'.format(sim.m[0]))
    text_base.append('\ndesign tankage length                      = {:.2f} m'.format(L))
    text_base.append('\ndesign mass flow rate                      = {:.2f} kg/s'.format(mdot))
    text_base.append('\ndesign airframe diameter                   = {:.2f} in.'.format(dia))
    text_base.append('\ndesign nozzle exit pressure                = {:.2f} kPa'.format(p_e))
    
    text_base.append('\n')
    text_base.append('\nCONSTRAINTS')
    text_base.append('\n-----------------------------')
    text_base.append('\nL/D ratio (c.f. < {})                      = {:.2f}'.format(cons_LD, ld_ratio(eng_sys_len, dia)))
    text_base.append('\nSommerfield criterion (c.f. pe/pa >= {})   = {:.1f}'.format(cons_S_crit, sim.S_crit))
    text_base.append("\nMax acceleration (c.f. < {})               = {:.2f} g's".format(cons_accel, max_g_force(sim.a)))
    text_base.append('\nTWR at lift off (c.f. > {})                = {:.2f}'.format(cons_TWR, sim.TWR))
    text_base.append('\naltitude at apogee (c.f. > {})             = {:.1f} km'.format(cons_alt/1000, sim.alt[-1]/1000))
    text_base.append('\nspeed when leaving launch rail (c.f. > {}) = {:.1f} m/s'.format(cons_ls, sim.launch_speed))
    text_base.append('\ndesign thrust (ground level) (c.f. < {})   = {:.1f} kN'.format(cons_thrust, sim.F[0]/1000))

    text_base.append('\n')
    text_base.append('\nADDITIONAL INFORMATION')
    text_base.append('\n-----------------------------')
    text_base.append('\naltitude at engine ignition                = {:.1f} m'.format(launch_site_alt))
    text_base.append('\nmission time at apogee                     = {:.1f} s'.format(sim.t[-1]))
    text_base.append('\nmission time at burnout                    = {:.1f} s'.format(sim.t[fdex-1]))
    text_base.append('\ndesign total propellant mass               = {:.3f} kg'.format(sim.m_prop[0]))
    text_base.append('\ndesign thrust (vacuum)                     = {:.1f} kN'.format(sim.F[fdex - 1]/1000))
    text_base.append('\ndesign burn time                           = {} s'.format(fdex*time_step))
    text_base.append('\ndesign expansion ratio                     = {:.1f}'.format(sim.ex))
    text_base.append('\ndesign throat area                         = {:.1f} in.^2'.format(sim.A_t/0.0254**2))
    text_base.append('\ndesign isp                                 = {:.1f} s'.format(sim.Ve/sim.g[0]))
    text_base.append('\ndesign chamber pressure                    = {:.1f} psi'.format(350))
    text_base.append('\ndesign total impulse                       = {:.1f} kN*s'.format(fdex*time_step*(sim.F[fdex - 1]/1000 + sim.F[0]/1000)/2))
    text_base.append('\ndesign dV                                  = {:.1f} km/s'.format(sim.dV1))
    text_base.append('\nestimated minimum required dV              = {:.1f} km/s'.format(sqrt(2*sim.g[0]*sim.alt[-1])/1000))
    return text_base

# this creates a nice set of plots of our trajectory data and saves it to rocket_farm
def rocket_plot(t, alt, v, a, F, q, Ma, m, p_a, D):
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import rc
    import pylab
    
    pylab.rcParams['figure.figsize'] = (10.0, 10.0)
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9) = plt.subplots(9, sharex=True)
    
    for n in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9):
        n.spines['top'].set_visible(False)
        n.spines['right'].set_visible(False)
        n.yaxis.set_ticks_position('left')
        n.xaxis.set_ticks_position('bottom')
        n.yaxis.labelpad = 20
        
    ax1.plot(t, alt/1000, 'k')
    ax1.set_ylabel("Altitude (km)")
    ax1.yaxis.major.locator.set_params(nbins=6)
    ax1.set_title('LV4 Trajectory')
    
    ax2.plot(t, v, 'k')
    ax2.yaxis.major.locator.set_params(nbins=6)
    ax2.set_ylabel("Velocity (m/s)")
    
    ax3.plot(t, a/g_0, 'k')
    ax3.yaxis.major.locator.set_params(nbins=10)
    ax3.set_ylabel("Acceleration/g_0")
    
    ax4.plot(t, F/1000, 'k')
    ax4.yaxis.major.locator.set_params(nbins=6)
    ax4.set_ylabel("Thrust (kN)")
    
    ax5.plot(t, q/1000, 'k')
    ax5.yaxis.major.locator.set_params(nbins=6)
    ax5.set_ylabel("Dynamic Pressure (kPa)")
    
    ax6.plot(t, Ma, 'k')
    ax6.yaxis.major.locator.set_params(nbins=6) 
    ax6.set_ylabel("Mach number")
    ax6.set_xlabel("t (s)")
    
    ax7.plot(t, np.array(m)*0.666*np.array(a), 'k')
    ax7.yaxis.major.locator.set_params(nbins=6) 
    ax7.set_ylabel("LOX Tank Axial Load")
    ax7.set_xlabel("t (s)")
    
    ax8.plot(t, D, 'k')
    ax8.yaxis.major.locator.set_params(nbins=6)
    ax8.set_ylabel("Drag (N)")
    
    ax9.plot(t, p_a/1000, 'k')
    ax9.yaxis.major.locator.set_params(nbins=6)
    ax9.set_ylabel("Air Pressure (Pa)")
    
    # we save the nice figures we make and then display them
    plt.savefig(openrkt.rkt_prefix +'psas_rocket_'+str(openrkt.get_index()-1)+'_traj.svg')
    plt.show()

# this creates some plots of the phase spaces of all our designs
def phase_plot(L, mdot, p_e):
    import pylab
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    
    ax1.plot(L, p_e)
    ax1.set_title('Design Vectors')
    ax1.yaxis.major.locator.set_params(nbins=6)
    ax1.set_xlabel("Length (m)")
    ax1.set_ylabel("Exit pressure (kPa)")
    
    ax2.plot(mdot, p_e)
    ax2.yaxis.major.locator.set_params(nbins=6)
    ax2.set_xlabel("Mass flow rate (kg/s)")
    ax2.set_ylabel("Exit pressure (kPa)")
    
    ax3.plot(L, mdot)
    ax3.yaxis.major.locator.set_params(nbins=6)
    ax3.set_ylabel("Mass flow rate (kg/s)")
    ax3.set_xlabel("Length (in)")
    
    # we display the first diagram of projected 2d phase portraits
    plt.show()
    
    fig2 = plt.figure()
    ax = fig2.add_subplot(111, projection='3d')
    
    ax.plot(L, mdot, p_e)
    ax.set_xlabel("Length (m)")
    ax.set_ylabel("Mass flow rate (kg/s)")
    ax.set_zlabel("Exit pressure (kPa)")
    
    # we display the interactive 3d phase portrait
    plt.show()
    # note, we're choosing to not automatically save these, but they can be saved from the interface
    

# Results, this is the big boi function
if __name__ == '__main__': # Testing  
    # feed initial design into iterative optimizer, get best design
    res = iterate(f, X0, iterations)
    print("Done!")
    
    # Rename the optimized output
    L = res[0]
    mdot = res[1]
    p_e = res[2]
    
    # get trajectory info from optimal design
    sim = trajectory.trajectory(L, mdot, dia, p_e, x_init=launch_site_alt, dt=time_step)  
    
    fdex = last_moment(sim.F) # index of engine burnout

    # get/print info about our trajectory and rocket
    res_text = print_results(res)
    for line in res_text:
        print(line)
    print('Engine system details in trajectory log!')
    print('\nMaking an OpenRocket rocket and corresponding engine!')
    # create an openrocket file with matching engine for our design (and print/save trajectory data)
    openrkt.make_engine(mdot, sim.m_prop[0], dia, sim.F[0:fdex], fdex*time_step, sim.Ve/sim.g[0], res_text)
    
    rocket_plot(sim.t, sim.alt, sim.v, sim.a, sim.F, sim.q, sim.Ma, sim.m, sim.p_a, sim.D) # draw pretty pictures of optimized trajectory
    
    # get us some nice plots of the phase space of design vectors
    designplot = [[],[],[]]
    for i in range(0, len(allvectors)):
        for j in range(0, len(designplot)):
            designplot[j].append(allvectors[i][j])
    phase_plot(designplot[0], designplot[1], designplot[2])
