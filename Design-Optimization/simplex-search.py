#Simplex search

#Erin Schmidt

#for non-linear programming problems ala Nelder and Mead(1965)
#**Note** use Python3 for this script
import math as m
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

def simplex_search(f, x_start, max_iter = 100, epsilon = 1E-6, gamma = 5, beta = 0.5):
	"""
	parameters of the function:
	f is the function to be optimized, return a scalar, and operate over a numpy array with the same dimensions of x_start
	x_start (numpy array): initial position
	epsilon is the termination criteria
	gamma is the contraction coefficient
	beta is the expansion coefficient
	"""
	
	#init arrays
	N = len(x_start)
	fnew = []
	xnew = []
	x = []
	
	#generate vertices of initial simplex
	a = .75
	x0 = (x_start)
	x1 = [x0 + [((N + 1)**0.5 + N - 1.)/(N + 1.)*a, 0.]]
	x2 = [x0 + [0., ((N + 1)**0.5 - 1.)/(N + 1.)*a]] 
	x3 = [x0 - [0., ((N + 1)**0.5 - 1.)/(N + 1.)*a]] 
	x = np.vstack((x1, x2, x3))
	#print(x)
	
	#simplex iteration
	while True:
		#find best, worst and 2nd worst points --> new center point
		f_run = np.array([f(x[0]), f(x[1]), f(x[2])]).tolist() #func. values at vertices
		#print(f_run)
		xw = x[f_run.index(sorted(f_run)[-1])]	#worst point
		xb = x[f_run.index(sorted(f_run)[0])] 	#best point
		xs = x[f_run.index(sorted(f_run)[-2])]	#2nd worst point
		xc = (xb + xs)/N #center point				  
		xr = 2*xc - xw #reflection point
		
		#check cases
		if f(xr) < f(xb): #expansion
			xnew = 2*xr - xc
			#xnew = (1 - gamma)*xc - gamma*xr
			#print('a', f(xr), f(xb)) #for debugging
		elif f(xr) > f(xw): #contraction 1
			xnew = (1 - beta)*xc + beta*xw
			#print('b', f(xr), f(xw))
		elif f(xs) < f(xr) and f(xr) < f(xw): #contraction 2
			xnew = (1 + beta)*xc - beta*xw
			#print('c', f(xs), f(xr), f(xw))
		else:
			xnew = xr
			
			#print('d', f(xr))
			#my xnews are total shit! They go in the wrong direction!
		
		#replace vertices
		"""
		if f(xnew) < f(xb):
		"""
		x[f_run.index(sorted(f_run)[-1])] = xnew
		#x[1] = xb
		#x[2] = xs
		fnew.append(f(xnew))
		print('Current optimum = ', fnew[-1])
		
		#break is any termination critera satisfied
		if len(fnew) == max_iter or term_check(xb, xc, xs, xnew, N) <= epsilon:
			return f(x[f_run.index(sorted(f_run)[0])]), x[f_run.index(sorted(f_run)[0])], len(fnew)
		"""if len(fnew)>1 and fnew[-1]-fnew[-2]<epsilon:
			return f(x[f_run.index(sorted(f_run)[0])]), x[f_run.index(sorted(f_run)[0])], len(fnew)"""

def term_check(xb, xc, xs, xnew, N): #the termination critera
	return m.sqrt(((f(xb) - f(xc))**2 + (f(xnew) - f(xc))**2 + (f(xs) - f(xc))**2)/(N + 1))

#testing
def f(z): #the objective function
	x=z
	return (1 - x[0])*(1 - x[0]) + (2 - x[1])*(2 - x[1])

#print results
(f, x, iter) = simplex_search(f, np.array([0,0]))
#print('\n')
print('f = ', f)
print('x = ', x)
print('iterations = ', iter)


#graphically verify the minimum
x = np.arange(-10,10,.1)
y = np.arange(-10,10,.1)
(x, y) = np.meshgrid(x, y)
z = (1 - x)**2 + (2 - y)**2
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(x, y, z, label='parametric curve', cmap=cm.jet, linewidth=0.2)
#ax.legend()
plt.show()
