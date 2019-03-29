#class simplex: # Nelder-Mead simplex search
import numpy as np
import trajectory
import liquid_motor

global cons_TWR, cons_S_crit, cons_accel, cons_LD, cons_alt, X0, m0

#CHANGE CONSTRAINTS HERE
cons_TWR = 3       # TWR constraint
cons_S_crit = 0.4#0.35 # Critical pressure ratio constraint
cons_accel = 15    # Max acceleration constraint
cons_LD = 15       # L/D ratio constraint
cons_alt = 100000  # Min altitude constraint
cons_F = 6

##CHANGE INITIAL DESIGN GUESS HERE
L = 1.47    # Tank length (m)
mdot = 2.393 #0.453592 * 0.9 * 5 # Propellant mass flow rate (kg/s)
dia = 12.204  # Rocket diameter (in)
p_e = 45.458  # Pressure (kPa)
X0 = np.array([L, mdot, dia, p_e])
m0 = trajectory.trajectory(X0[0], X0[1], X0[2], X0[3])[-4][0] # Initial guess GLOW

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
            (alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop) = sim.trajectory(xb[0], xb[1], xb[2], xb[3])
            return f(x[f_run.index(sorted(f_run)[0])], rp), x[f_run.index(sorted(f_run)[0])], len(fb)
        
def term_check(N, rp, f_run, fxc): # Termination critera
    M = [0]*(N + 1)
    for i in range(0, N + 1):
        if i == f_run.index(sorted(f_run)[-1]): # Avoid worst point
            M[i] = 0
        else:
            M[i] = (f_run[i] - fxc)**2
    #return m.sqrt(((f(xb) - f(xc))**2 + (f(xnew) - f(xc))**2 + (f(xs) - f(xc))**2)/(N + 1))
    return m.sqrt(sum(M)/N)
        
# Pseudo-objective function
def f(x, rp=25):
    #x = np.array(x)
    L = x[0]   # Rocket length (m)
    mdot = x[1] # Propellant mass flow rate (kg/s)
    dia = x[2] # Rocket diameter (in)
    p_e = x[3]  # Pressure (kPa)
    (alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop) = trajectory.trajectory(L, mdot, dia, p_e)
    obj_func = m[0]/m0 + rp*(max(0, F[0]/(cons_F*1000) - 1)**2 + max(0, (L+2)/(dia*0.0254*cons_LD) - 1)**2 + \
               max(0, -TWR/cons_TWR + 1)**2 + max(0, -S_crit/cons_S_crit + 1)**2 + \
                             max(0, -alt[-1]/cons_alt + 1)**2 + max(0, max(abs(a))/(cons_accel*9.81) - 1)**2)   
    #print "obj: ", obj_func, " alt: ",alt[-1]/1000, " F: ",F[0]/1000,"  L, mdot, dia, p_e: ", x
    return obj_func

# Results
if __name__ == '__main__': # Testing
    from math import sqrt, pi, exp, log, cos
    import math as m
    
    """max_iter = 200
    rp = 50
    gamma = 6
    beta = .5
    a = 5
    (f, x, it) = search(f, np.array(X0), max_iter, gamma, beta, rp, a)
    """
    
    from scipy.optimize import minimize
    #res = optimize.basinhopping(f, X0)    
    #res = minimize(f, X0, method='nelder-mead', options={'adaptive':True})
    res = minimize(f, X0, method='Powell', options={'disp': True})
    
    (alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop) = trajectory.trajectory(res.x[0], res.x[1], res.x[2], res.x[3]) #, 500)   
    print('\n')
    
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import rc
    import pylab
    #%config InlineBackend.figure_formats=['svg']
    #%matplotlib inline
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    
    # Rename the optimized output
    L = res.x[0]
    mdot = res.x[1]
    dia = res.x[2]
    p_e = res.x[3]
    
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
    #plt.savefig('traj.svg')
    #plt.show()
    
    
    np.set_printoptions(precision=3)
    print('\n')
    print('OPTIMIZED DESIGN VECTOR')
    print('-----------------------------')
    print('x_optimized                                = ', res.x)
    print('x_initial_guess                            = ', X0)
    print('design tankage length                      = {0:.2f} m'.format(res.x[0]))
    print('design mass flow rate                      = {0:.2f} kg/s'.format(res.x[1]))
    print('design airframe diameter                   = {0:.2f} in.'.format(res.x[2]))
    print('design nozzle exit pressure                = {0:.2f} kPa'.format(res.x[3]))
    print('iterations                                 =', res.nit)
    print('design GLOW                                = {0:.1f} kg'.format(m[0]))
    print('x0 GLOW                                    = {0:.1f} kg'.format( \
          trajectory.trajectory(X0[0], X0[1], X0[2], X0[3])[-4][0]))

    print('\n')
    print('CONSTRAINTS')
    print('-----------------------------')
    print('L/D ratio (c.f. < {})                      = {:.2f}'.format(cons_LD, (L+2)/(dia*0.0254)))
    print('Sommerfield criterion (c.f. pe/pa >= {})   = {:.1f}'.format(cons_S_crit, S_crit))
    print("Max acceleration (c.f. < {})               = {:.2f} g's".format(cons_accel, max(abs(a))/9.81))
    print('TWR at lift off (c.f. > {})                = {:.2f}'.format(cons_TWR, TWR))
    print('altitude at apogee                         = {0:.1f} km'.format(alt[-1]/1000))

    print('\n')
    print('ADDITIONAL INFORMATION')
    print('-----------------------------')
    print('mission time at apogee                     = {0:.1f} s'.format(t[-1]))
    print('design total propellant mass               = {0:.3f} kg'.format(m_prop[0]))
    print('design thrust (sea level)                  = {0:.1f} kN'.format(F[0]/1000))
    j = 0
    for thing in F:
        if thing == 0:
            fdex = j
            break
        j += 1
    print('design thrust (vacuum)                     = {0:.1f} kN'.format(F[fdex - 1]/1000))
    print('design burn time                           = {} s'.format(fdex/10))
    print('design expansion ratio                     = {0:.1f}'.format(ex))
    print('design throat area                         = {0:.1f} in.^2'.format(A_t/0.0254**2))
    print('design isp                                 = {0:.1f} s'.format(Ve/9.81))
    print('design chamber pressure                    = {0:.1f} psi'.format(350))
    print('design total impulse                       = {0:.1f} kN*s'.format(j/10*(F[fdex - 1]/1000 + F[0]/1000)/2))
    print('design dV                                  = {0:.1f} km/s'.format(dV1))
    print('estimated minimum required dV              = {0:.1f} km/s'.format(sqrt(2*9.81*alt[-1])/1000))
    print('\n Making an OpenRocket engine!')
    liquid_motor.make_engine(mdot, m_prop[0], dia, F[fdex - 1], fdex/10, Ve/9.81)
