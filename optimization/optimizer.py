# could improve:
# constraints
# initial design vectors
# rocket evaluator
# minimizer
# overall redundancy reduction

# Nelder-Mead simplex search, could also use Powell's method

import numpy as np
import trajectory
import openrocket_interface as openrkt
from math import sqrt, pi, exp, log, cos
import math as m # is this sanitary? is this needed? 'm' is a variable name too...
from scipy.optimize import minimize
global cons_TWR, cons_S_crit, cons_accel, cons_LD, cons_alt, X0, m0 # i'm sorry for global vars...
global allvectors, cons_mass

# SIMULATION AND OPTIMIZATION PARAMETERS
time_step = 0.1 # change time-step for trajectory simulation
iterations = 3 # number of escalating iterations, degenerate after ~8
launch_site_alt = 1401 # m, altitude of launch site above sea level

##CHANGE INITIAL DESIGN GUESS HERE
# these are kinda janky guesses, but they let you see designs get sucked into what seems to be an attractor
L = 1.7    # Total tank lengths (m)
mdot = 2.37 # Propellant mass flow rate (kg/s)
dia = 12.32  # Rocket diameter (in)
p_e = 45.56  # Exit Pressure (kPa)

#CHANGE CONSTRAINTS HERE
###CHECK ME
cons_mass = trajectory.trajectory(L, mdot, dia, p_e, x_init=launch_site_alt, dt=time_step)[-5][0] # GLOW constraint taken from initial design GLOV
cons_totimp = 400   # total impulse constraint
cons_ls = 20.       # min launch speed from 60' tower constraint, m/s
cons_TWR = 1.1       # TWR constraint
cons_S_crit = 0.35 # Critical pressure ratio constraint
cons_accel = 15.    # Max acceleration constraint
cons_LD = 20.       # L/D ratio constraint
cons_alt = 100000 + launch_site_alt # Min altitude constraint
cons_thrust = 6.    # max average thrust
max_dia = 14.    # maximum diameter

# physical constants
nose_l = 1.25 # m
ers_l = 6 * 0.0254 # m (converted from in)
rcs_l = 6 * 0.0254 # m (converted from in)
av_l = 18 * 0.0254 # m (converted from in)
n2_l = 18 * 0.0254 # m (converted from in)

shirt_l = sum([nose_l, ers_l, rcs_l, av_l, n2_l]) # m, total of above. (what about gaps??)

g_0 = openrkt.g_0 # gravity

X0 = np.array([L, mdot, dia, p_e]) #numpy arrays are nicer
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

def at_most(var, cons):
    return (var - cons)
    #return max(0, var/cons - 1)**2
    
def at_least(var, cons):
    return -(at_most(var, cons))
    #return max(0, 1 - var/cons)**2
    
def obj_func(glow):
    return glow/cons_mass
    
def merit_func(thrust, rho):
    return (rho/2) * (thrust - cons_thrust)**2

def pen_func(ls, alt, S_crit, LD, gees, TWR, mu):
    print(LD)
    print(cons_LD)
    return mu * sum([log(-at_least(ls, cons_ls)), log(-at_least(alt, cons_alt)), \
        log(-at_least(S_crit, cons_S_crit)), max(0, LD/cons_LD - 1)**2, \
        log(-at_most(gees, cons_accel)), log(-at_least(TWR, cons_TWR))])

# basic exterior penalty function
# evaluate a rocket design and trajectory
# (a/b - 1): a < b
# (1 - a/b): a > b
def eval_rkt(ls, thrust, LD, TWR, S_crit, alt, gees, dia, punisher=25):
    # multiply all constraint violations, degenerate around 10^9
    combined = punisher*(\
    
    # max launch speed
    1*max(0, 1 - ls/cons_ls)**2 + \
    
    # min engine thrust
    1*max(0, thrust/cons_thrust - 1)**2 + \
    
    # min tank length / airframe diameter ratio
    1*max(0, LD/cons_LD - 1)**2 + \
    
    # max thrust / weight ratio
    1*max(0, 1 - TWR/cons_TWR)**2 + \
    
    # max Sommerfield constant constraint
    1*max(0, 1 - S_crit/cons_S_crit)**2 + \
    
    # max altitude
    1*max(0, 1 - alt/cons_alt)**2 + \
    
    # min top g's
    1*max(0, gees/cons_accel - 1)**2 + \
    1*max(0, dia/max_dia - 1)**2) 
    
    return combined

# Pseudo-objective function
# x is array of design parameters, n is degree of penalty function
def f(x, n):
    global allvectors
    L = x[0]   # Tank length (m)
    mdot = x[1] # Propellant mass flow rate (kg/s)
    dia = x[2] # Rocket diameter (in)
    p_e = x[3]  # Pressure (kPa)
    
    # get trajectory data from x
    (alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop, ls) = trajectory.trajectory(L, mdot, dia, p_e, x_init=launch_site_alt, dt=time_step)
    
    obj = m[0]/cons_mass # use this line to minimize mass
    
    # or use 3 of the next 4 lines of code to maximize total impulse
    fdex = last_moment(F) # index of engine burnout
    #total_impulse = (fdex * time_step) * (F[fdex - 1] + F[0]) / 2 / 1000 # burn-time times average thrust in kN*s
    
    #obj_func = cons_totimp / total_impulse # this is one way to maximize
    #obj_func = -total_impulse / cons_totimp # this is another way to maximize
    
    # then, calculate penalization from trajectory based on initial thrust
    r, l_o, l_f = openrkt.split_tanks(m[0], dia)
    pen = eval_rkt(ls, F[fdex-1]/1000, ld_ratio(openrkt.system_length(l_o, l_f), dia), TWR, S_crit, alt[-1], max_g_force(a), dia, 10**n)
    
    # add objective and penalty functions
    #sum_func = obj_func + pen_func
    #obj = obj_func(m[0]) + merit_func(max(F)/1000, 10)
    #pen = pen_func(ls, alt[-1], S_crit, ld_ratio(L, dia), max_g_force(a), TWR, 10)
    sum_func = obj + pen
    
    # print blocks are sanity checks so i'm not staring at a blank screen, commented out because it's unnecessary'
    #print(L, mdot, dia, p_e, alt[-1])
    print(obj, ' + ', pen, ' = ', sum_func)
    
    allvectors.append(x) # maintains a list of every successive design
    
    return sum_func

# we want to iterate our optimizer to get nice results  
# n = number of iterations, it gets weird for n > 8, I keep it around 2-6  
def iterate(f, x_0, n):
    x = x_0 # initial design vector
    
    for i in range(n):
        print("Iteration " + str(i+1) + ":")
        #res = minimize(f, x, args=(i+1), method='Powell', options={'disp': True}) # this is a minimizer for brute-force monkeys
        res = minimize(f, x, args=(i+1), method='nelder-mead', options={'disp': True, 'adaptive':True}) # this minimizer uses simplex method
        x = res.x # feed optimal design vec into next iteration
        #cons_mass = trajectory.trajectory(x[0], x[1], x[2], x[3], x_init=launch_site_alt, dt=time_step)[-5][0] # update mass constraint
        print("Propellant tube length (m): "+str(x[0]))
        print("Mass flow rate (kg/s): "+str(x[1]))
        print("Airframe diameter (in): "+str(x[2]))
        print("Exit pressure (kPa): "+str(x[3]))
    return x

# this creates a list of strings for relevant data of trajectory
def print_results(res):
    text_base = [] # list of lines of strings
    
    # Rename the optimized output
    L = res[0]
    mdot = res[1]
    dia = res[2]
    p_e = res[3]
    
    (alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop, ls) = trajectory.trajectory(L, mdot, dia, p_e, x_init=launch_site_alt, dt=time_step) 
    fdex = last_moment(F) # index of engine burnout
    
    np.set_printoptions(precision=3) # this line may be deprecated, i copy-pasted most of this section
    
    text_base.append('\nOPTIMIZED DESIGN VECTOR')
    text_base.append('\n-----------------------------')
    text_base.append('\nx_initial_guess                            = ' + ', '.join([str(X0[0]), str(X0[1]), str(X0[2]), str(X0[3])])) # kind of a kludge
    text_base.append('\ninitial guess GLOW                                    = {:.1f} kg'.format( \
          trajectory.trajectory(X0[0], X0[1], X0[2], X0[3], x_init=launch_site_alt, dt=time_step)[-5][0]))
    text_base.append('\nx_optimized                                = ' + ', '.join([str(L), str(mdot), str(dia), str(p_e)])) # kind of a kludge
    text_base.append('\ndesign GLOW                                = {:.1f} kg'.format(m[0]))
    text_base.append('\ndesign tankage length                      = {:.2f} m'.format(L))
    text_base.append('\ndesign mass flow rate                      = {:.2f} kg/s'.format(mdot))
    text_base.append('\ndesign airframe diameter                   = {:.2f} in.'.format(dia))
    text_base.append('\ndesign nozzle exit pressure                = {:.2f} kPa'.format(p_e))
    
    text_base.append('\n')
    text_base.append('\nCONSTRAINTS')
    text_base.append('\n-----------------------------')
    text_base.append('\nL/D ratio (c.f. < {})                      = {:.2f}'.format(cons_LD, ld_ratio(L, dia)))
    text_base.append('\nSommerfield criterion (c.f. pe/pa >= {})   = {:.1f}'.format(cons_S_crit, S_crit))
    text_base.append("\nMax acceleration (c.f. < {})               = {:.2f} g's".format(cons_accel, max_g_force(a)))
    text_base.append('\nTWR at lift off (c.f. > {})                = {:.2f}'.format(cons_TWR, TWR))
    text_base.append('\naltitude at apogee (c.f. > {})             = {:.1f} km'.format(cons_alt/1000, alt[-1]/1000))
    text_base.append('\nspeed when leaving launch rail (c.f. > {}) = {:.1f} km/s'.format(cons_ls, ls))
    text_base.append('\ndesign thrust (vacuum) (c.f. < {})         = {:.1f} kN'.format(cons_thrust, F[fdex - 1]/1000))

    text_base.append('\n')
    text_base.append('\nADDITIONAL INFORMATION')
    text_base.append('\n-----------------------------')
    text_base.append('\naltitude at engine ignition                = {:.1f} m'.format(launch_site_alt))
    text_base.append('\nmission time at apogee                     = {:.1f} s'.format(t[-1]))
    text_base.append('\ndesign total propellant mass               = {:.3f} kg'.format(m_prop[0]))
    text_base.append('\ndesign thrust (ground level)                  = {:.1f} kN'.format(F[0]/1000))
    text_base.append('\ndesign burn time                           = {} s'.format(fdex*time_step))
    text_base.append('\ndesign expansion ratio                     = {:.1f}'.format(ex))
    text_base.append('\ndesign throat area                         = {:.1f} in.^2'.format(A_t/0.0254**2))
    text_base.append('\ndesign isp                                 = {:.1f} s'.format(Ve/g_0))
    text_base.append('\ndesign chamber pressure                    = {:.1f} psi'.format(350))
    text_base.append('\ndesign total impulse                       = {:.1f} kN*s'.format(fdex*time_step*(F[fdex - 1]/1000 + F[0]/1000)/2))
    text_base.append('\ndesign dV                                  = {:.1f} km/s'.format(dV1))
    text_base.append('\nestimated minimum required dV              = {:.1f} km/s'.format(sqrt(2*g_0*alt[-1])/1000))
    
    return text_base

# this creates a nice set of plots of our trajectory data and saves it to rocket_farm
def rocket_plot(t, alt, v, a, F, q, Ma, m):
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import rc
    import pylab
    #%config InlineBackend.figure_formats=['svg']
    #%matplotlib inline
    #rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    #rc('text', usetex=True)
    
    pylab.rcParams['figure.figsize'] = (10.0, 10.0)
    fig, (ax1, ax2, ax3, ax4, ax6, ax7, ax8) = plt.subplots(7, sharex=True)
    
    for n in (ax1, ax2, ax3, ax4, ax6, ax7, ax8):
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
    
    ax6.plot(t, q/1000, 'k')
    ax6.yaxis.major.locator.set_params(nbins=6)
    ax6.set_ylabel("Dynamic Pressure (kPa)")
    
    ax7.plot(t, Ma, 'k')
    ax7.yaxis.major.locator.set_params(nbins=6) 
    ax7.set_ylabel("Mach number")
    ax7.set_xlabel("t (s)")
    
    ax8.plot(t, np.array(m)*0.666*np.array(a), 'k')
    ax8.yaxis.major.locator.set_params(nbins=6) 
    ax8.set_ylabel("LOX Tank Axial Load")
    ax8.set_xlabel("t (s)")
    
    # we save the nice figures we make and then display them
    plt.savefig(openrkt.rkt_prefix +'psas_rocket_'+str(openrkt.get_index()-1)+'_traj.svg')
    plt.show()

# this creates some plots of the phase spaces of all our designs
def phase_plot(L, mdot, D, p_e):
    import pylab
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    
    ax1.plot(L, D)
    ax1.set_title('Design Vectors')
    ax1.yaxis.major.locator.set_params(nbins=6)
    ax1.set_xlabel("Length (m)")
    ax1.set_ylabel("Diameter (in)")
    
    ax2.plot(mdot, p_e)
    ax2.yaxis.major.locator.set_params(nbins=6)
    ax2.set_xlabel("Mass flow rate (kg/s)")
    ax2.set_ylabel("Exit pressure (kPa)")
    
    ax3.plot(mdot, L)
    ax3.yaxis.major.locator.set_params(nbins=6)
    ax3.set_xlabel("Mass flow rate (kg/s)")
    ax3.set_ylabel("Length (m)")
    
    # we display the first diagram of 2d phase portraits
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
    dia = res[2]
    p_e = res[3]
    
    # get trajectory info from optimal design
    (alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop, ls) = trajectory.trajectory(L, mdot, dia, p_e, x_init=launch_site_alt, dt=time_step)  
    
    fdex = last_moment(F) # index of engine burnout

    # get/print info about our trajectory and rocket
    res_text = print_results(res)
    for line in res_text:
        print(line)
 
    print('\nMaking an OpenRocket rocket and corresponding engine!')
    # create an openrocket file with matching engine for our design (and print/save trajectory data)
    openrkt.make_engine(mdot, m_prop[0], dia, F[0:fdex], fdex*time_step, Ve/g_0, res_text)
    
    rocket_plot(t, alt, v, a, F, q, Ma, m) # draw pretty pictures
    
    # the ugly code below is to get us some nice plots of the phase space of design vectors
    y0, y1, y2, y3 = [], [], [], []
    for i in range(0, len(allvectors)):
        y0.append(allvectors[i][0])
        y1.append(allvectors[i][1])
        y2.append(allvectors[i][2])
        y3.append(allvectors[i][3])
    phase_plot(y0, y1, y2, y3)
