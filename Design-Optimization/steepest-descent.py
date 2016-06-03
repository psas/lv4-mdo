#Steepest descent alorithem
#-Erin Schmidt

#
import numpy as np
import math as m

class eq_int:
	def func(alpha, x, norm_del_f, rp, count, n=0): 		#the objective function  	
		count += 1
		x[0] += alpha*norm_del_f[0]
		x[1] += alpha*norm_del_f[1]
		f = st_dec.func(x, rp, n)[0]
		return f, count
	
	def mini(au, al, x, norm_del_f, rp, count): 	#evaluates f at the minimum (or optimum) stationary point
		alpha = (au + al)*0.5
		(f, count) = eq_int.func(alpha, x, norm_del_f, rp, count)
		return f, alpha, count
    
	def search(al, x, norm_del_f, rp, delta=0.01, epsilon=1E-4, count=0):
		(f, count) = eq_int.func(al, x, norm_del_f, rp, count)
		fl = f 					#function value at lower bound
		#delta = 0.01 			#step-size
		#au = 0.15 				#alpha upper bound

		while True:
			aa = delta
			(f, count) = eq_int.func(aa, x, norm_del_f, rp, count)
			fa = f
			if fa > fl:
				delta = delta * 0.1
			else:
				break

		while True:
			au = aa + delta
			(f, count) = eq_int.func(au, x, norm_del_f, rp, count)
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
					(f, count) = eq_int.func(au, x, norm_del_f, rp, count)
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
				(f, alpha, count) = eq_int.mini(au, al, x, norm_del_f, rp, count)
				return f, alpha, count
       
class st_dec:
    def func(x, rp, n=0):
        g = max([0, 170 - x[0]])
        f = (x[0] - 100)**2 + (x[1] - 50)**2 + rp*g**2 #pseudo-objective value at x
        del_f = [2*(x[0] - 100) + rp*(-2*g), 2*(x[1] -50)] #gradient value at x
        
        #Test function w/o constraint
        #f = (x[0] -2)**4 + (x[0] - 2*x[1])**2
        #del_f = [4*(x[0] - 2)**3 + 2*(x[0] -2*x[1]), -4*(x[0] -2*x[1])]
        del_f = np.array(del_f)
        f_val = (x[0] - 100)**2 + (x[1] - 50)**2 #actual function value?
        n += 1
        return f, del_f, n, f_val
    
    def steepest(x, rp):
        n = 0 #iteration counter
        alpha = .01 #initial step size
        (f, del_f, n, f_val) = st_dec.func(x, rp, n)
        while True:
            #f_old = f
            alpha_old = alpha
            norm_del_f = -del_f/m.sqrt(del_f[0]**2+del_f[1]**2) #normalize the gradient vector
            alpha = eq_int.search(alpha, x, norm_del_f, rp)[1]
            #print(alpha)
            x[0] += alpha*norm_del_f[0] #next step x-values
            x[1] += alpha*norm_del_f[1]  
            (f, del_f, n, f_val) = st_dec.func(x, rp, n)
            if abs((alpha - alpha_old)/alpha) < 1E-10: #convergence criteria
                return f, n, x, alpha, f_val
            if n > 100000:
                return f, n, x, alpha, f_val

#starting design        
x = [0,0]
rp = 100
(f, n, x, alpha, f_val) = st_dec.steepest(x, rp)
print(f_val, x, n)