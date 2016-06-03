#Equal Interval Search
#Erin Schmidt

#Adapted, with significant modification, from Arora et al.'s APOLLO implementation found in "Introduction to Optimum Design" 1st Ed. (1989).
from __future__ import division
import numpy as np
import math as m

class eq_int:
	def func(alpha, count=0): 		#the objective function
		count += 1
		#f = x**2 - 10*x +34
		x = [0,3]
		del_f = [4*(x[0]-2)**3 + 2*(x[0] -2*x[1]), -4*(x[0] - 2*x[1])]
		del_f = np.array(del_f)
		norm_del_f = -del_f/m.sqrt(del_f[0]**2+del_f[1]**2)
		x[0] += norm_del_f[0]*alpha
		x[1] += norm_del_f[1]*alpha
		f = (x[0] - 2)**4 + (x[0] - 2*x[1])**2
      #good, looks like it works!
		return f, count
	
	def mini(au, al, count): 	#evaluates f at the minimum (or optimum) stationary point
		alpha = (au + al)*0.5
		(f, count) = eq_int.func(alpha, count)
		return f, alpha, count
    
	def search(delta=0.01, epsilon=1E-3, count=0):
		al = 0.01  			    #should this be zero? alpha lower bound
		(f, count) = eq_int.func(al, count)
		fl = f 					#function value at lower bound
		#delta = 0.01 			#step-size
		#au = 0.15 				#alpha upper bound

		while True:
			aa = delta
			(f, count) = eq_int.func(aa, count)
			fa = f
			if fa > fl:
				delta = delta * 0.1
			else:
				break

		while True:
			au = aa + delta
			(f, count) = eq_int.func(au, count)
			fu = f
			if fa > fu:
				al = aa
				aa = au
				fl = fa
				fa = fu
			else:
				break

		while True:
			if (au - al) > epsilon: #compares interval size to convergence criteria
				delta = delta * 0.1
				aa = al 			#intermediate alpha
				fa = fl 			#intermediate alpha function value
				while True:
					au = aa + delta
					(f, count) = eq_int.func(au, count)
					fu = f
					if fa > fu: 
						al = aa
						aa = au
						fl = fa
						fa = fu
						continue
					else:
						break
				continue
			else:
				(f, alpha, count) = eq_int.mini(au, al, count)
				return f, alpha, count

(f, alpha, count) = eq_int.search()
print('The minimum is at {:.4f}'.format(alpha))
print('The function value at the minimum = {:.4f}'.format(f))
print('Total number of function calls = {}'.format(count))
	
