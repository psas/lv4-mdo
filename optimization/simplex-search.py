class optim: # Nelder-Mead simplex search
#-Erin Schmidt

    import math as m
    import numpy as np
    
    def simplex_search(f, x_start, max_iter = 100, epsilon = 1E-6, gamma = 5, beta = 0.5):
        
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
        fnew = []        # Empty function matrix
        xnew = []        # Empty re-write for design variables
        x    = []        # Empty x matrix
        C    = [0]*N        # Empty center point matrix
        
        # Generate vertices of initial simplex
        a = .75          # Step Size Alpha
        x0 = (x_start)   # x0 Value for x Matrix
        x1 = [x0 + [((N + 1)**0.5 + N - 1.)/(N + 1.)*a, 0., 0., 0.]]
        x2 = [x0 + [0., ((N + 1)**0.5 - 1.)/(N + 1.)*a, 0., 0.]]
        x3 = [x0 + [0., 0., ((N + 1)**0.5 - 1.)/(N + 1.)*a, 0.]]
        x4 = [x0 + [0., 0., 0., ((N + 1)**0.5 - 1.)/(N + 1.)*a]]
        x = optim.np.vstack((x0, x1, x2, x3, x4))
        #print(x)
    
        # Simplex iteration
        while True:
            # Find best, worst, 2nd worst, and new center point
            f_run = optim.np.array([optim.f(x[0]), optim.f(x[1]), optim.f(x[2]), optim.f(x[3]), optim.f(x[4])]).tolist() # Func. values at vertices
            #print(f_run)
            xw = x[f_run.index(sorted(f_run)[-1])] # Worst point
            xb = x[f_run.index(sorted(f_run)[0])]  # Best point
            xs = x[f_run.index(sorted(f_run)[-2])] # 2nd worst point        
            # xc = (xb + xs)/N                     # Center point
            for i in range(0, N):
                if optim.f(x[i]) == optim.f(xw):
                    C[i] = 0
                else:
                    C[i] = x[i]
            xc = (sum(C)/N)                        # Center point
            xr = 2*xc - xw                         # Reflection point
            
            # Check cases
            if optim.f(xr) < optim.f(xb):                      # Expansion
                xnew = (1 + gamma)*xc - gamma*xr
            elif optim.f(xr) > optim.f(xw):                    # Contraction 1
                xnew = (1 - beta)*xc + beta*xw
            elif optim.f(xs) < optim.f(xr) and optim.f(xr) < optim.f(xw):  # Contraction 2
                xnew = (1 + beta)*xc - beta*xw
            else:
                xnew = xr
            
            # Replace Vertices
            x[f_run.index(sorted(f_run)[-1])] = xnew
            # x[1] = xb
            # x[2] = xs
            fnew.append(f(xnew))
            # print('Current optimum = ', fnew[-1])
            
            # Break if any termination critera is satisfied
            if len(fnew) == max_iter or optim.term_check(x, xc, xw, N) <= epsilon:
                return optim.f(x[f_run.index(sorted(f_run)[0])]), x[f_run.index(sorted(f_run)[0])], len(fnew)
            
    def term_check(x, xc, xw, N): # Termination critera
        M = [0]*N
        for i in range(0, N):
            if optim.f(x[i]) == optim.f(xw): # Avoid worst point
                M[i] = 0
            else:
                M[i] = (optim.f(x[i]) - optim.f(xc))**2
        # return m.sqrt(((f(xb) - f(xc))**2 + (f(xnew) - f(xc))**2 + (f(xs) - f(xc))**2 + ()**2)/(N + 1))
        return m.sqrt(sum(M)/(N+1))
        
    # Pseudo-objective function
    def f(x): 
        rp = 100
        L = x[0]   # Rocket length (m)
        mdot = x[1] # Propellant mass flow rate (kg/s)
        dia = x[2] # Rocket diameter (in)
        p_e = x[3]  # Pressure (kPa)
        (L, dia, alt, v, a, t, F, D, Ma, rho, p_a, T_a, TWR, ex, Ve, A_t, dV1, m, S_crit, q, m_prop) = sim.trajectory(L, mdot, dia, p_e)
        obj_func = m[0] + rp*(max(0, (L+2)/(dia*0.0254) - 15)**2 + max(0, -TWR + 2)**2 + max(0, -S_crit + 0.35)**2 + max(0, -alt[-1] +100000)**2 + max(0, max(a)/9.81 - 15)**2)
        return obj_func

# Results
from math import sqrt, pi, exp, log, cos
import numpy as np
import math as m

(f, x, iter) = optim.simplex_search(optim.f, np.array([1.5, 0.453592 * 0.9 * 4, 8, 50]))
print('\n')
print('f = ', f)
print('x = ', x)
print('iterations = ', iter)

