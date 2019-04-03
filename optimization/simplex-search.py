# could improve:
# constraints
# initial design vectors
# rocket evaluator
# minimizer
# figure plotting/saving
# overall redundancy reduction

#class simplex: # Nelder-Mead simplex search

import numpy as np
import trajectory
import openrocket_interface as openrkt
from math import sqrt, pi, exp, log, cos
import math as m # is this sanitary? is this needed? 'm' is a variable name too...
from scipy.optimize import minimize
global cons_TWR, cons_S_crit, cons_accel, cons_LD, cons_alt, X0, m0 # i'm sorry for global vars...

#CHANGE CONSTRAINTS HERE, these have been changed from the MDO, need sanity check from ME's
time_step = .1 # change time-step for trajectorys
cons_mass = 400     # GLOW constraint
cons_TWR = 1.5       # TWR constraint
cons_S_crit = 0.35  #0.35 # Critical pressure ratio constraint
cons_accel = 15    # Max acceleration constraint
cons_LD = 18       # L/D ratio constraint
cons_alt = 115000  # Min altitude constraint
cons_thrust = 5.5         # max average thrust
cons_ls = 25.       # min launch speed, m/s

##CHANGE INITIAL DESIGN GUESS HERE
L = 1.979   # Tank length (m)
mdot = 2.318 #0.453592 * 0.9 * 5 # Propellant mass flow rate (kg/s)
dia = 12.392  # Rocket diameter (in)
p_e = 54.202  # Pressure (kPa)
#X0 = np.array([L, mdot, dia, p_e]) #idk if this needs to be a np.array
X0 = [L, mdot, dia, p_e]
m0 = trajectory.trajectory(X0[0], X0[1], X0[2], X0[3], dt=time_step)[-5][0] # Initial guess GLOW

# this function is deprecated pls ignore
def search(f, x_start, max_iter = 100, gamma = 5, beta = 0.5, rp=100, a=10, epsilon = 1E-6):
    
    """
    parameters of the function:
    f is the function to be optimized
    x_start (numpy array) is the initial simplex vertices
    epsilon is the termination criteria
    gamma is the contraction coefficient
    beta is the expansion coefficient
    """
    
    # Init Arrays
    N = len(x_start) # Amount of design variables
    fb = []        # Empty function matrix
    xnew = []        # Empty re-write for design variables
    x    = []        # Empty x matrix
    C    = [[0]*N]*(N+1)        # Empty center point matrix #####CHANGED
    
    # Generate vertices of initial simplex
    x0 = (x_start)   # x0 Value for x Matrix
    x1 = [x0 + [((N + 1)**0.5 + N - 1.)/(N + 1.)*a, 0., 0., 0.]]
    x2 = [x0 + [0., ((N + 1)**0.5 - 1.)/(N + 1.)*a, 0., 0.]]
    x3 = [x0 + [0., 0., ((N + 1)**0.5 - 1.)/(N + 1.)*a, 0.]]
    x4 = [x0 + [0., 0., 0., ((N + 1)**0.5 - 1.)/(N + 1.)*a]]
    x = np.vstack((x0, x1, x2, x3, x4))

    # Simplex iteration
    while True:
        # Find best, worst, 2nd worst, and new center point
        f_run = np.array([f(x[0], rp), f(x[1], rp), f(x[2], rp), f(x[3], rp), f(x[4], rp)]).tolist() # Func. values at vertices
        xw = x[f_run.index(sorted(f_run)[-1])] # Worst point
        xb = x[f_run.index(sorted(f_run)[0])]  # Best point
        xs = x[f_run.index(sorted(f_run)[-2])] # 2nd worst point        
        # xc = (xb + xs)/N                     # Center point
        for i in range(0, N+1):
            if i == f_run.index(sorted(f_run)[-1]):
                C[i] = [0,0,0,0]
            else:
                C[i] = x[i].tolist()
        xc = sum(np.array(C))/(N) # Center point
        xr = 2*xc - xw # Reflection point
        fxr = f(xr, rp)
        fxc = f(xc, rp)
        
        # Check cases
        if fxr < f_run[f_run.index(sorted(f_run)[0])]: #f(xr, rp) < f(xb, rp): # Expansion
            xnew = (1 + gamma)*xc - gamma*xr
        elif fxr > f_run[f_run.index(sorted(f_run)[-1])]: #f(xr, rp) > f(xw, rp): # Contraction 1
            xnew = (1 - beta)*xc + beta*xw
        elif f_run[f_run.index(sorted(f_run)[-2])] < fxr and fxr < f_run[f_run.index(sorted(f_run)[-1])]: #f(xs, rp) < f(xr, rp) and f(xr, rp) < f(xw, rp): # Contraction 2
            xnew = (1 + beta)*xc - beta*xw
        else:
            xnew = xr
        
        # Replace Vertices
        x[f_run.index(sorted(f_run)[-1])] = xnew
        #x[f_run.index(sorted(f_run)[1])] = xb
        #x[f_run.index(sorted(f_run)[2])] = xs
        fb.append(f(xb, rp))
        print('Current optimum = ', fb[-1])
        
        # Break if any termination critera is satisfied
        if len(fb) == max_iter: # or term_check(x, xc, xw, N, rp, f_run) <= epsilon:
            (alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop, ls) = sim.trajectory(xb[0], xb[1], xb[2], xb[3])
            return f(x[f_run.index(sorted(f_run)[0])], rp), x[f_run.index(sorted(f_run)[0])], len(fb)

# this function is deprecated pls ignore        
def term_check(N, rp, f_run, fxc): # Termination critera
    M = [0]*(N + 1)
    for i in range(0, N + 1):
        if i == f_run.index(sorted(f_run)[-1]): # Avoid worst point
            M[i] = 0
        else:
            M[i] = (f_run[i] - fxc)**2
    #return m.sqrt(((f(xb) - f(xc))**2 + (f(xnew) - f(xc))**2 + (f(xs) - f(xc))**2)/(N + 1))
    return m.sqrt(sum(M)/N)

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

# basic exterior penalty function
# evaluate a rocket design and trajectory
# (a/b - 1): a < b
# (1 - a/b): a > b
def eval_rkt(ls, thrust, LD, TWR, S_crit, alt, gees, punisher=25):
    # multiply all constraint violations
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
    1*max(0, gees/cons_accel - 1)**2) 
    return combined

# Pseudo-objective function
# x is array of design parameters, n is degree of penalty function
def f(x, n):
    x = np.array(x)
    L = x[0]   # Rocket length (m)
    mdot = x[1] # Propellant mass flow rate (kg/s)
    dia = x[2] # Rocket diameter (in)
    p_e = x[3]  # Pressure (kPa)
    # get trajectory data from x
    (alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop, ls) = trajectory.trajectory(L, mdot, dia, p_e, dt=time_step)
    #obj_func = m[0]/cons_mass # use this to minimize mass
    #obj_func = 400/(Ve/9.8066) # or use this to maximize specific impulse
    
    # or use this block to maximize total impulse
    fdex = last_moment(F)
    impulse = fdex/10*(F[fdex - 1]/1000 + F[0]/1000)/2
    obj_func = (400/impulse)
    
    pen_func = eval_rkt(ls, F[0]/1000, (L+2)/(dia*0.0254), TWR, S_crit, alt[-1], max(abs(a))/9.81, 10^n)
    sum_func = obj_func + pen_func
    # print blocks are sanity checks so i'm not staring at a blank screen
    print(L, mdot, dia, p_e, alt[-1])
    print(obj_func, ' + ', pen_func, ' = ', sum_func)
    return sum_func

# we want to iterate our optimizer to get nice results, n is number of iterations
# it gets weird for n > 8, I keep it around 3-7    
def iterate(f, x_0):
    n=1
    x = x_0
    for i in range(n):
        print(i, ": \n")
        #res = minimize(f, x, args=(i+1), method='Powell', options={'disp': True}) # this is a minimizer for monkeys
        res = minimize(f, x, args=(i+1), method='nelder-mead', options={'disp': True, 'adaptive':True}) # this minimizer may or may not suck
        x = res.x
    return res

# this creates a list of strings for relevant data of trajectory
def print_results(res):
    text_base = []
    # Rename the optimized output
    L = res.x[0]
    mdot = res.x[1]
    dia = res.x[2]
    p_e = res.x[3]
    
    (alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop, ls) = trajectory.trajectory(L, mdot, dia, p_e, dt=time_step) 
    fdex = last_moment(F)
    np.set_printoptions(precision=3) # this line may be deprecated, i copy-pasted most of this section
    
    text_base.append('\nOPTIMIZED DESIGN VECTOR')
    text_base.append('\n-----------------------------')
    text_base.append('\nx_optimized                                = ' + ', '.join([str(L), str(mdot), str(dia), str(p_e)])) # kind of a kludge
    text_base.append('\nx_initial_guess                            = ' + ', '.join([str(X0[0]), str(X0[1]), str(X0[2]), str(X0[3])])) # kind of a kludge
    text_base.append('\ndesign tankage length                      = {:.2f} m'.format(L))
    text_base.append('\ndesign mass flow rate                      = {:.2f} kg/s'.format(mdot))
    text_base.append('\ndesign airframe diameter                   = {:.2f} in.'.format(dia))
    text_base.append('\ndesign nozzle exit pressure                = {:.2f} kPa'.format(p_e))
    text_base.append('\niterations (of last iteration)             = ' + str(res.nit))
    text_base.append('\ndesign GLOW                                = {:.1f} kg'.format(m[0]))
    text_base.append('\nx0 GLOW                                    = {:.1f} kg'.format( \
          trajectory.trajectory(X0[0], X0[1], X0[2], X0[3], dt=time_step)[-5][0]))

    text_base.append('\n')
    text_base.append('\nCONSTRAINTS')
    text_base.append('\n-----------------------------')
    text_base.append('\nL/D ratio (c.f. < {})                      = {:.2f}'.format(cons_LD, (L+2)/(dia*0.0254)))
    text_base.append('\nSommerfield criterion (c.f. pe/pa >= {})   = {:.1f}'.format(cons_S_crit, S_crit))
    text_base.append("\nMax acceleration (c.f. < {})               = {:.2f} g's".format(cons_accel, max(abs(a))/9.81))
    text_base.append('\nTWR at lift off (c.f. > {})                = {:.2f}'.format(cons_TWR, TWR))
    text_base.append('\naltitude at apogee (c.f. > {})             = {:.1f} km'.format(cons_alt, alt[-1]/1000))
    text_base.append('\nspeed when leaving launch rail (c.f. > {}) = {:.1f} km/s'.format(cons_ls, ls))
    text_base.append('\ndesign thrust (vacuum) (c.f. < {})         = {:.1f} kN'.format(cons_thrust, F[fdex - 1]/1000))

    text_base.append('\n')
    text_base.append('\nADDITIONAL INFORMATION')
    text_base.append('\n-----------------------------')
    text_base.append('\nmission time at apogee                     = {:.1f} s'.format(t[-1]))
    text_base.append('\ndesign total propellant mass               = {:.3f} kg'.format(m_prop[0]))
    text_base.append('\ndesign thrust (sea level)                  = {:.1f} kN'.format(F[0]/1000))
    
    text_base.append('\ndesign burn time                           = {} s'.format(fdex/10))
    text_base.append('\ndesign expansion ratio                     = {:.1f}'.format(ex))
    text_base.append('\ndesign throat area                         = {:.1f} in.^2'.format(A_t/0.0254**2))
    text_base.append('\ndesign isp                                 = {:.1f} s'.format(Ve/9.81))
    text_base.append('\ndesign chamber pressure                    = {:.1f} psi'.format(350))
    text_base.append('\ndesign total impulse                       = {:.1f} kN*s'.format(fdex/10*(F[fdex - 1]/1000 + F[0]/1000)/2))
    text_base.append('\ndesign dV                                  = {:.1f} km/s'.format(dV1))
    text_base.append('\nestimated minimum required dV              = {:.1f} km/s'.format(sqrt(2*9.81*alt[-1])/1000))
    return text_base

# Results, this is the big boi function
if __name__ == '__main__': # Testing  
    #this block is deprecated
    """max_iter = 200
    rp = 50
    gamma = 6
    beta = .5
    a = 5
    (f, x, it) = search(f, np.array(X0), max_iter, gamma, beta, rp, a)
    """
    #res = optimize.basinhopping(f, X0) # idk why this is here
    
    # feed initial design into iterative optimizer, get best design
    res = iterate(f, X0)
    
    # Rename the optimized output
    L = res.x[0]
    mdot = res.x[1]
    dia = res.x[2]
    p_e = res.x[3]
    
    # get trajectory info from optimal design
    (alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop, ls) = trajectory.trajectory(L, mdot, dia, p_e, dt=time_step)  
    fdex = last_moment(F)
    print('\n')
    print("params: ", res)
    
    # get info about our trajectory and rocket
    res_text = print_results(res)
    for line in res_text:
        print(line)
 
    print('\n Making an OpenRocket engine!')
    # create an openrocket file for our design (and print trajectory data!)
    openrkt.make_engine(mdot, m_prop[0], dia, F[0], fdex/10, Ve/9.81, res_text)
    
    # functional code ends here

    # this matplotlib nonsense is from the og script, nothing below here matters to me
    # i just ignore it all, since it crashes if showing/saving are uncommented
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import rc
    import pylab
    #%config InlineBackend.figure_formats=['svg']
    #%matplotlib inline
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    
    pylab.rcParams['figure.figsize'] = (10.0, 10.0)
    f, (ax1, ax2, ax3, ax4, ax6, ax7, ax8) = plt.subplots(7, sharex=True)
    for n in (ax1, ax2, ax3, ax4, ax6, ax7, ax8):
        n.spines['top'].set_visible(False)
        n.spines['right'].set_visible(False)
        n.yaxis.set_ticks_position('left')
        n.xaxis.set_ticks_position('bottom')
        n.yaxis.labelpad = 20
        
    ax1.plot(t, alt/1000, 'k')
    ax1.set_ylabel("Altitude (km)")
    ax1.yaxis.major.locator.set_params(nbins=6)
    #ax1.set_title('LV4 Trajectory')
    ax2.plot(t, v, 'k')
    ax2.yaxis.major.locator.set_params(nbins=6)
    ax2.set_ylabel("Velocity (m/s)")
    ax3.plot(t, a/9.81, 'k')
    ax3.yaxis.major.locator.set_params(nbins=10)
    ax3.set_ylabel("Acceleration/g0")
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
    #ax8.yaxis.major.locator.set_params(nbins=6) 
    ax8.set_ylabel("LOX Tank Axial Load")
    ax8.set_xlabel("t (s)")
    #plt.savefig('traj.svg') # this crashes
    #plt.show() # this crashes
