#First put problem in standard form and construct cannonical matrix
#     ^^ program will NOT do this.

#Note: this is coded in Python3 *not* MATLAB

 #initial feasible solution is given by initial tableu construction

from __future__ import division

def printTableu(tableu):
    print('----------------------------------------------')  
    
    for row in tableu:
        newrow = [ '%2.2f' % elem for elem in row ]
        print(newrow)
    print('----------------------------------------------')
    return


def pivoter(tableu, row, col): #where row and col are the index values of min b[i]/a[i][j]
    j = 0
    pivot = tableu[row][col]
    for x in tableu[row]: #normalize entire tableu by the pivot value
        tableu[row][j] = tableu[row][j] / pivot
        j += 1
        
    i = 0
    for xi in tableu: #Gauss-Jordan elimination
        if i != row: #ignore the pivot row
            ratio = xi[col]
            j = 0
            for xij in xi: #subtract by columns
                xij -= ratio * tableu[row][j]
                tableu[i][j] = xij
                j += 1
                
        i += 1 #move across rows
        
    return tableu

def simplex(tableu):
     THETA_INFINITE = -1
     optimal = False
     unbounded  = False
     n = len(tableu[3])
     m = len(tableu) - 1
     while ((not optimal) and (not unbounded)):
         min = 0.0
         pivotCol = j = 0
         while(j < (n-m)): #find min of C[j]
              cj = tableu[3][j]
              if (cj < min) and (j > 0):
                  min = cj
                  #pivotCol = j
              j += 1   
         if min == 0.0: #if not C[j] < 0 then break the loop
             optimal = True
             continue
         
         pivotRow = i = 0
         minTheta = THETA_INFINITE
         for xi in tableu:
              if (i > 0):
                   xij = xi[pivotCol]
                   if xij > 0: #to avoid infinite and negative numbers
                       theta = (xi[0] / xij) #test criteria for pivot column -> min b[i]/a[i][j]
                       if (theta < minTheta) or (minTheta == THETA_INFINITE):
                           minTheta = theta
                           pivotRow = i
                           
              i += 1
         if minTheta == THETA_INFINITE:
             unbounded = True
             continue
         
         tableu = pivoter(tableu, pivotRow, pivotCol)
     
     print('\n Unbounded = {}'.format(unbounded))
     print('Optimal = {} \n'.format(optimal))
     print("Final tableu")
     printTableu(tableu)
     return tableu

f  = [ 0.0, -1.0, -1.0, -2.0, 0.0,  0.0]
x1 = [ 8.0, 2.0,  1.0,  2.0, -1.0,  0.0]
x2 = [ 2.0, 1.0,  1.0,  1.0, 0.0, 1.0]
x3 = [ 1.0, -1.0,  1.0,  2.0, 0.0,  0.0]


tableu = []
tableu.append(x1)
tableu.append(x2)
tableu.append(x3)
tableu.append(f)

print("Initial tableu")
printTableu(tableu)

tableu = simplex(tableu)