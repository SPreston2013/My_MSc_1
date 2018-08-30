import matplotlib.pyplot as plt
import os
from scipy.misc import derivative
import math 
import numpy as np


Nlist = np.linspace(50, 56.8, 10) #Change in e-folds

#Create empty arrays 
Phi_g_array = []
f_r_store = []
f_O_store = []
f_ns_store = []

nlist  = np.linspace(0.865,2,60) #Powers of n

          

#Newton to solve for phi_gw
def Newton(f, dfdx, x, eps):
    f_value = f(x, N, n)
    iteration_counter = 0
    while abs(f_value) > eps and iteration_counter < 100:
        try:
            x = x - float(f_value)/dfdx(x, n)
        except ZeroDivisionError:
            print ("Error! - derivative zero for x = ", x)
            sys.exit(1)    

        f_value = f(x, N, n)
        iteration_counter += 1

  
    if abs(f_value) > eps:
        iteration_counter = -1
    return x, iteration_counter

#Equation for r
def f_r(N, n) :
    return 16.*n/(4.*N + n)

#Equation for Omega_gw
def f_O(N, n, g) : 
    return 4.36e-15*r*(np.sqrt(2.*N*n + 0.5*n**2)/g)**(-n) 

#Equation for n_s
def f_ns(N, n) :
    return 1. - (2.0*(n + 2.0)/(4.*N + n))


#Quadratic to be solved by newton to find phi_gw
def f(x, N, n):
    return 1/(2.0*n)*x**2 + (n/2.0)*math.log(x) - n/4. - (n/2.0)*math.log((np.sqrt(2.*N*n + (n**2/2.0)))) - 37.09829

#The derivative of the above quadratic 
def dfdx(x, n):
    return x/n + 0.5*n/x                                                                        

#loop for values of n and N                                                                    
for j in range(len(nlist)):
    n = nlist[j]

    for i in range(len(Nlist)):
            N = Nlist[i]
            g, no_iterations  = Newton(f, dfdx, x=1000, eps=1.0e-6) #Call newton for phi_gw

            r = f_r(N, n)         
            O = f_O(N, n, g)
            ns = f_ns(N, n) 
            f_r_store.append(r)
            f_O_store.append(O)
            f_ns_store.append(ns)
        

    if no_iterations > 0:
        #print ("Number of function calls: %d" % (1 + 2*no_iterations))
        print ("The value of phi_gw is: %f" % (g))
    else:
        print ("Solution not found!")



 #   plt.plot(f_r_store, f_O_store, color = 'red')  #this plots Omega_gw against r
    plt.plot(f_ns_store, f_O_store, color = 'blue')  #this plots Omega_gw against ns

    f_r_store=[]
    f_O_store=[] 
    f_ns_store=[] 



    
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
#plt.xlabel(r'$r$',fontsize=22)
#plt.axvline(x=0.07)
plt.xlabel(r'$n_s$',fontsize=18)
plt.axvline(x=0.9591)
plt.axvline(x=0.9749) 
plt.ylabel(r'$\Omega_{gw}h^2$',fontsize=18)


plt.show()
    
