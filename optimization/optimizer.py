# could improve:
# minimizer
# actual analysis: error, sensitivity, stability

# Nelder-Mead simplex search
# 6 iterations with timestep 0.25 took my laptop 6 minutes to complete
# need to: error analysis on trajectory to find ideal time-step
#          somehow analyze optimizer to find proper number of iterations
#          which might require writing our own optimizer that won't attempt
#          undefined math operations outside the barrier

import numpy as np
from math import sqrt, pi, exp, log, cos
from scipy.optimize import minimize, differential_evolution

import trajectory
import openrocket_interface as openrkt

# i'm sorry for global vars... i'll sanitize namespace soonish
global cons_TWR, cons_S_crit, cons_accel, cons_LD, cons_alt, X0, m0, cons_mass 
global allvectors, dbz
dbz = 0 # divisions by zero (from crossing constraint boundary)

# SIMULATION AND OPTIMIZATION PARAMETERS
time_step = 0.25       # change time-step for trajectory simulation, need to find best step size
iterations = 6         # number of escalating iterations, degenerates after 6!
launch_site_alt = 1401 # m, altitude of launch site above sea level

##CHANGE INITIAL DESIGN GUESS HERE
# be sure that you start with a feasible design, otherwise the problem will be ill-conditioned
L = 1.5    # Total tank lengths (m)
mdot = 2.2 # Propellant mass flow rate (kg/s)
dia = 12.  # Rocket diameter (in)
p_e = 47.  # Exit Pressure (kPa)

#CHANGE CONSTRAINTS HERE
cons_mass = 200.                         # GLOW constraint, kg
cons_ls = 22.                            # min launch speed from 60' tower constraint, m/s
cons_TWR = 2.                            # TWR constraint
cons_S_crit = 0.35                       # Critical pressure ratio constraint
cons_accel = 15.                         # Max acceleration constraint, g's
cons_LD = 18.                            # L/D ratio constraint
cons_alt = 100000. + launch_site_alt     # Min altitude constraint, km
cons_thrust = 6.                         # max ground-level thrust, kN
cons_ceiling = 150000. + launch_site_alt # base-11 maximum apogee requirement, km

shirt_l = trajectory.shirt_l  # non-engine subsystems lengths
g_n = openrkt.g_n             # standard gravity
X0 = np.array([L, mdot, p_e]) # numpy arrays are nicer
allvectors = []               # array for all design vecs, global variable

# calculates length-diameter ratio
def ld_ratio(pants_l, dia):
    total_l = shirt_l + pants_l # add shirt and pants lengths for total rocket length
    D = dia * 0.0254            # converts in to m for total airframe diameter
    return total_l / D
    
# calculates top g-force
def max_g_force(a):
    return max(abs(a)) / g_n
    
# all of our comparisons are ratios instead of subtractions because
# it's normalized, instead of dependant on magnitudes of constraints

# minimize this, **2 makes it well behaved w.r.t. when var=cons
def objective(var, cons):
    return (var/cons)**2 / 2

# **2 because i like it more than abs(), but that also works
def exact(var, cons):
    return (var/cons - 1)**2 / 2

# this is your basic exterior penalty, either punishes for unfeasibility or is inactive
def exterior(var, cons, good_if_less_than=False):
    if good_if_less_than:
        return max(0, var/cons - 1)**2 / 2
    else:
        return max(0, -(var/cons - 1))**2 / 2

# this restricts our objective function to the strictly feasible region
# make rockets great again, build that wall, etc, watch out for division by zero
# i like the logarithmic version more, but interior point method is viable too
# i'm not sure if our try/except block effects the theoretical convergence result
# theoretically, the barriers should just work but the simplex method is semi-blind
def barrier(var, cons, int_point=False, good_if_less_than=True):
    global dbz
    def interior(g): # in case we don't like logarithms
        return 1/g
    try:
        #print('log: '+str(-log((var/cons - 1)))+ ', 1/x: '+str(-interior(-(var/cons - 1)))) # debug
        if not int_point:
            if good_if_less_than:
                return -log(-(var/cons - 1))
            else:
                return -log((var/cons - 1))
        elif int_point:
            if good_if_less_than:
                return -interior(var/cons - 1)
            else:
                return -interior(-(var/cons - 1))
    except:
        dbz += 1 # keep track of divisions by zero, side effect
        return 10**10

# this manages all our constraints
# penalty parameters: mu -> 0 and rho -> infinity 
def penalty(ls, F, LD, TWR, S_crit, alt, max_g, mu, rho):
    b = [barrier(alt, cons_alt, int_point=False, good_if_less_than=False),
         barrier(alt, cons_ceiling, int_point=False, good_if_less_than=True)]
    eq = [exact(F, cons_thrust)]
    ext = [exterior(ls, cons_ls, good_if_less_than=False),
           exterior(LD, cons_LD, good_if_less_than=True),
           exterior(TWR, cons_TWR, good_if_less_than=False),
           exterior(S_crit, cons_S_crit, good_if_less_than=False),
           exterior(max_g, cons_accel, good_if_less_than=True)]
    #print("barrier: " + str(mu*sum(b))+", exact: "+ str(rho/2 *sum(eq))+", exterior: "+ str(rho*sum(ext))) # debug
    return mu*sum(b) + rho*(sum(eq) + sum(ext))

# Pseudo-objective function
# x is array of design parameters, n is sequence index of penalty and barrier functions
# print blocks are sanity checks so i'm not staring at a blank screen and can see what various tweaks actually do
def f(x, n=4):
    global allvectors
    L    = x[0]  # Tank length (m)
    mdot = x[1]  # Propellant mass flow rate (kg/s)
    p_e  = x[2]  # Pressure (kPa)
    
    # get trajectory data from x
    sim = trajectory.trajectory(L, mdot, dia, p_e, x_init=launch_site_alt, dt=time_step)
    r, l_o, l_f = openrkt.split_tanks(sim.m_prop[0], dia)
    
    obj_func = objective(sim.m[0], cons_mass) # minimize GLOW
    # then, calculate penalization from trajectory based on initial thrust
    pen_func = penalty(sim.launch_speed,
                  sim.F[0]/1000,
                  ld_ratio(openrkt.system_length(l_o, l_f), dia),
                  sim.TWR,
                  sim.S_crit,
                  sim.alt[-1],
                  max_g_force(sim.a),
                  (0.80) / (2**n),   # initial mu and
                  (1.25) * (2**n)) #                rho are selected for nice behavior
    # add objective and penalty functions
    merit_func = obj_func + pen_func
    
    # DEBUG BLOCK
    #print("obj: "+str(obj))
    #print("Total (n): "+str(sum_func) + ' ('+ str(n) +')')
    #print("L "+ str(L)+ ", D " + str(dia)+ ", mdot " + str(mdot)+ ", p_e "+  str(p_e)+ ", alt "+  str(sim.alt[-1]))
    #print("obj "+ str(obj)+ ' + pen '+ str(pen)+ ' = '+ str(sum_func))
    #print('')
    allvectors.append(x) # maintains a list of every design, side effect
    return merit_func

# we want to iterate our optimizer for theoretical convergence reasons (given some assumptions)
# n = number of sequential iterations, one day this will terminate itself automatically 
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
        print("Exit pressure (kPa): "+str(x[2]))
        print("Altitude (km): "+str(alt[-1]))
        print('')
        dbz=0 # I only care about divisions by zero in each individual iteration, side effect
    return x

# this creates a list of strings for relevant data of trajectory
def print_results(res):
    text_base = [] # list of lines of strings
    
    # Rename the optimized output for convenience
    L    = res[0]
    mdot = res[1]
    p_e  = res[2]
    
    sim = trajectory.trajectory(L, mdot, dia, p_e, x_init=launch_site_alt, dt=time_step) 
    r, l_o, l_f = openrkt.split_tanks(sim.m_prop[0], dia)
    eng_sys_len = openrkt.system_length(l_o, l_f)
    
    np.set_printoptions(precision=3) # this line may be deprecated, i copy-pasted most of this section
    
    text_base.append('\nOPTIMIZED DESIGN VECTOR')
    text_base.append('\n-----------------------------')
    text_base.append('\nx_initial_guess                            = ' + ', '.join([str(X0[0]), str(X0[1]), str(X0[2])]))
    text_base.append('\ninitial guess GLOW                                    = {:.1f} kg'.format( \
          trajectory.trajectory(X0[0], X0[1], dia, X0[2], x_init=launch_site_alt, dt=time_step).m[0]))
    text_base.append('\nx_optimized                                = ' + ', '.join([str(L), str(mdot), str(p_e)]))
    text_base.append('\ndesign tankage length                      = {:.3f} m'.format(L))
    text_base.append('\ndesign mass flow rate                      = {:.3f} kg/s'.format(mdot))
    text_base.append('\ndesign nozzle exit pressure                = {:.3f} kPa'.format(p_e))
    text_base.append('\ndesign airframe diameter                   = {:.3f} in.'.format(dia))
    text_base.append('\ndesign GLOW                                = {:.3f} kg'.format(sim.m[0]))
    
    text_base.append('\n')
    text_base.append('\nCONSTRAINTS')
    text_base.append('\n-----------------------------')
    text_base.append('\nL/D ratio (c.f. < {})                      = {:.3f}'.format(cons_LD, ld_ratio(eng_sys_len, dia)))
    text_base.append('\nSommerfield criterion (c.f. pe/pa >= {})   = {:.3f}'.format(cons_S_crit, sim.S_crit))
    text_base.append("\nMax acceleration (c.f. < {})               = {:.3f} g's".format(cons_accel, max_g_force(sim.a)))
    text_base.append('\nTWR at lift off (c.f. > {})                = {:.3f}'.format(cons_TWR, sim.TWR))
    text_base.append('\naltitude at apogee (c.f. > {})             = {:.3f} km'.format(cons_alt/1000, sim.alt[-1]/1000))
    text_base.append('\nspeed when leaving launch rail (c.f. > {}) = {:.3f} m/s'.format(cons_ls, sim.launch_speed))
    text_base.append('\ndesign thrust (ground level) (c.f. < {})   = {:.3f} kN'.format(cons_thrust, sim.F[0]/1000))

    text_base.append('\n')
    text_base.append('\nADDITIONAL INFORMATION')
    text_base.append('\n-----------------------------')
    text_base.append('\naltitude at engine ignition                = {:.1f} m'.format(launch_site_alt))
    text_base.append('\nmission time at apogee                     = {:.3f} s'.format(sim.t[-1]))
    text_base.append('\nmission time at burnout                    = {:.3f} s'.format(sim.t[sim.F_index-1]))
    text_base.append('\ndesign total propellant mass               = {:.3f} kg'.format(sim.m_prop[0]))
    text_base.append('\ndesign thrust (vacuum)                     = {:.2f} kN'.format(sim.F[sim.F_index - 1]/1000))
    text_base.append('\ndesign burn time                           = {} s'.format(sim.F_index*time_step))
    text_base.append('\ndesign expansion ratio                     = {:.3f}'.format(sim.ex))
    text_base.append('\ndesign throat area                         = {:.3f} in.^2'.format(sim.A_t/0.0254**2))
    text_base.append('\ndesign isp                                 = {:.3f} s'.format(sim.Ve/g_n))
    text_base.append('\ndesign chamber pressure                    = {:.3f} psi'.format(350))
    text_base.append('\ndesign total impulse                       = {:.3f} kN*s'.format(sim.F_index*time_step*(sim.F[sim.F_index - 1]/1000 + sim.F[0]/1000)/2))
    text_base.append('\ndesign dV                                  = {:.3f} km/s'.format(sim.dV1))
    text_base.append('\nestimated minimum required dV              = {:.3f} km/s'.format(sqrt(2*g_n*sim.alt[-1])/1000))
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
    
    ax3.plot(t, a/g_n, 'k')
    ax3.yaxis.major.locator.set_params(nbins=10)
    ax3.set_ylabel("Acceleration/g_n")
    
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

# this creates some plots of the phase spaces of all our designs, doesn't save them
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
    
    # this block is for probing design space within bounds
    # these parameters were an experiment, idk how genetic algorithms work at all
    #res = differential_evolution(f, [(1.2, 1.9), (1.5, 3.), (30, 50)], \
                #strategy='best1exp', popsize=50, mutation=(.6,1.8), recombination=.1, polish=True,workers=-1, disp=True)
    #res=res.x
    
    print("Done!")
    
    # Rename the optimized output for convenience
    L    = res[0]
    mdot = res[1]
    p_e  = res[2]
    
    # get trajectory info from optimal design
    sim = trajectory.trajectory(L, mdot, dia, p_e, x_init=launch_site_alt, dt=time_step)  
    
    # get/print info about our trajectory and rocket
    res_text = print_results(res)
    for line in res_text:
        print(line)
    print('Engine system details in trajectory log!')
    print('\nMaking an OpenRocket rocket and corresponding engine!')
    
    # create an openrocket file with matching engine for our design (and print/save trajectory data)
    openrkt.make_engine(mdot, sim.m_prop[0], dia, sim.F[0:sim.F_index], sim.F_index*time_step, sim.Ve/g_n, res_text)
    
    # draw pretty pictures of optimized trajectory
    rocket_plot(sim.t, sim.alt, sim.v, sim.a, sim.F, sim.q, sim.Ma, sim.m, sim.p_a, sim.D)
    
    # get us some nice plots of the phase space of design vectors
    designplot = [[],[],[]]
    for i in range(0, len(allvectors)):
        for j in range(0, len(designplot)):
            designplot[j].append(allvectors[i][j])
    phase_plot(designplot[0], designplot[1], designplot[2])
