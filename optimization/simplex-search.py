#class simplex: # Nelder-Mead simplex search
import numpy as np

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
def f(x, rp=50):
    #x = np.array(x)
    L = x[0]   # Rocket length (m)
    mdot = x[1] # Propellant mass flow rate (kg/s)
    dia = x[2] # Rocket diameter (in)
    p_e = x[3]  # Pressure (kPa)
    (alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop) = sim.trajectory(L, mdot, dia, p_e)
    obj_func = m[0] + rp*(max(0, (L+2)/(dia*0.0254) - 15)**2 + max(0, -TWR + 2)**2 + max(0, -S_crit + 0.35)**2 + max(0, -alt[-1] + 100000)**2 + max(0, max(abs(a))/9.81 - 15)**2)   
    #obj_func = m[0] #+ rp*(max(0, -alt[-1] + 100000)**2)    
    return obj_func

# Results
if __name__ == '__main__': # Testing
    import numpy as np
    from trajectory import sim
    from math import sqrt, pi, exp, log, cos
    import math as m
    
    #X0 = np.array([1, 0.453592 * 0.9 * 4, 12, 50])
    X0 = np.array([2, 0.453592 * 0.9 * 6, 8, 50])
    """max_iter = 200
    rp = 50
    gamma = 6
    beta = .5
    a = 5
    (f, x, it) = search(f, np.array(X0), max_iter, gamma, beta, rp, a)
    """
    
    from scipy.optimize import minimize
    #res = optimize.basinhopping(f, X0)    
    res = minimize(f, X0, method='nelder-mead')    
    
    
    (alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop) = sim.trajectory(res.x[0], res.x[1], res.x[2], res.x[3])   
    print('\n')
    
    import matplotlib
    import matplotlib.pyplot as plt
    import pylab
    #%config InlineBackend.figure_formats=['svg']
    #%matplotlib inline
    
    L = x[0]
    mdot = x[1]
    dia = x[2]
    p_e = x[3]
    
    pylab.rcParams['figure.figsize'] = (10.0, 10.0)
    f, (ax1, ax2, ax3, ax4, ax6, ax7) = plt.subplots(6, sharex=True)
    #plt.xlim(0, 1.8)
    ax1.plot(t, alt/1000)
    ax1.set_ylabel("Altitude (km)")
    ax1.yaxis.major.locator.set_params(nbins=6)
    ax1.set_title('LV4 Trajectory')
    ax2.plot(t, v)
    ax2.yaxis.major.locator.set_params(nbins=6)
    ax2.set_ylabel("Velocity (m/s)")
    ax3.plot(t, a/9.81)
    ax3.yaxis.major.locator.set_params(nbins=10)
    ax3.set_ylabel("Acceleration/g0")
    ax4.plot(t, F/1000)
    ax4.yaxis.major.locator.set_params(nbins=6)
    ax4.set_ylabel("Thrust (kN)")
    ax6.plot(t, q/1000)
    ax6.yaxis.major.locator.set_params(nbins=6)
    ax6.set_ylabel("Dynamic Pressure (kPa)")
    ax7.plot(t, Ma)
    ax7.yaxis.major.locator.set_params(nbins=6) 
    ax7.set_ylabel("Mach number")
    ax7.set_xlabel("t (s)")
    plt.show()
    
    print('\n')
    print('MINIMIZED VALUE')
    print('-----------------------------')
    print('x = ', res.x)
    print('x0 = ', X0)
    print('iterations = ', it)
    print('design GLOW = {0:.1f} kg'.format(m[0]))
    print('x0 GLOW = ', sim.trajectory(X0[0], X0[1], X0[2], X0[3])[-4][0])
    
    print('\n')
    print('CONSTRAINTS')
    print('-----------------------------')
    print('L/D ratio (check < 15) = {0:.2f}'.format((L+2)/(dia*0.0254)))
    print('Sommerfield criterion (check pe/pa >= 0.3) = {0:.1f}'.format(S_crit))
    print('Max acceleration (check < 15) = {0:.2f}'.format(max(abs(a))/9.81))
    print('TWR at lift off (check TWR > 2) = {0:.2f}'.format(TWR))
    print('altitude at apogee = {0:.1f} km'.format(alt[-1]/1000))
    
    print('\n')
    print('ADDITIONAL INFORMATION')
    print('-----------------------------')
    print('mission time at apogee = {0:.1f} s'.format(t[-1]))
    print('design total propellant mass = {0:.3f}'.format(m_prop[0]))
    print('design thrust (sea level) = {0:.1f} kN'.format(F[0]/1000))
    j = 0
    for i in F:
        if i == 0:
            fdex = j
            break
        j += 1
    print('design thrust (vacuum) = {0:.1f} kN'.format(F[fdex - 1]/1000))
    print('design burn time = {} s'.format(fdex))
    print('design expansion ratio = {0:.1f}'.format(ex))
    print('design throat area = {0:.1f} in^2'.format(A_t/0.0254**2))
    print('design isp = {0:.1f} s'.format(Ve/9.81))
    #print('design dV = {} km/s c.f. required potential energy est = {} km/s'.format(dV1, sqrt(2*9.81*alt[-1])/1000))
    