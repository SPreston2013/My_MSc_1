import getdist.plots as gplot
import os
from pylab import *
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

G=gplot.getSinglePlotter(chain_dir=r'/Users/scottpreston/Desktop/base')
roots = ['base_r_plikHM_TTTEEE_lowTEB']
pairs = [[u'ns', u'r']]
G.plots_2d(roots, param_pairs=pairs, filled=True, shaded=False)
ylim(ymax=0.25)


N = np.linspace(50.0, 60.0, 10) 
n = 2/3 
n_1 = 2
n_2 = 4/3
n_3 = 3
n_4 = 1
n_5 = 1.125

def f_r0(N, n): # phi 2/3
    for j in range(len(N)):
        return (16*n)/(4*N + n)

def ns_N0(N, n): 
    for k in range(len(N)):
        return 1 -2*(n + 2)/(4*N + n)

def f_r1(N): #R^2 inflation 
    for j in range(len(N)):
        return 12/N**2

def ns_N1(N): 
    for k in range(len(N)):
        return 1 - 2/N 

def f_r2(N, n_1): #Massive potential 
    for j in range(len(N)):
        return (16*n_1)/(4*N + n_1)

def ns_N2(N, n_1): 
    for k in range(len(N)):
        return 1 -2*(n_1 + 2)/(4*N + n_1)

def f_r3(N, n_2): # phi 4/3
    for j in range(len(N)):
        return (16*n_2)/(4*N + n_2)

def ns_N3(N, n_2): 
    for k in range(len(N)):
        return 1 -2*(n_2 + 2)/(4*N + n_2)

def f_r4(N, n_3): # phi 3
    for j in range(len(N)):
        return (16*n_3)/(4*N + n_3)

def ns_N4(N, n_3): 
    for k in range(len(N)):
        return 1 -2*(n_3 + 2)/(4*N + n_3)

def f_r5(N, n_4): # phi 
    for j in range(len(N)):
        return (16*n_4)/(4*N + n_4)

def ns_N5(N, n_4): 
    for k in range(len(N)):
        return 1 -2*(n_4 + 2)/(4*N + n_4)

def f_r6(N, n_5): #Random potential 
    for j in range(len(N)):
        return (16*n_5)/(4*N + n_5)

def ns_N6(N, n_5): 
    for k in range(len(N)):
        return 1 -2*(n_5 + 2)/(4*N + n_5)

####################################################    
Nlist = np.linspace(50.0, 60.0, 10)
mulist = np.logspace(1.05,3, 120)            
N_1=0        
mu=0

phi_end_array = []   
f_r_store  = []       
f_ns_store = []       

def newton(f, x, dfdx, epsilon=1.0E-7, N_1=100, store=False):
    f_value = f(x)
    n=0
    if store : info= [(x,f_value)]
    while abs(f_value) > epsilon and n<=N_1:
        dfdx_value = float(dfdx(x))
        if abs(dfdx_value) < 1E-14:
            raise ValueError("Newton: f'(%g)=%g" % (x,dfdx_value))
        x = x -f_value/dfdx_value
        n += 1
        f_value = f(x)
        if store: info.append((x, f_value))
    if store:
        return x,info
    else:
        return x,n,f_value

#####################################
def f_r(phistar) :
        return 128.0*(phistar**3/((1.0-phistar**4)*mu))**2

def f_ns(phistar) : 
        return -24.0*phistar**2/(mu**2*(1.0-phistar**4))-0.375*f_r(phistar)+1.0

#####################################
def g(x):                   
    c = np.sqrt(8.0)/mu
    return  x**4 +c*x**3 -1.0   # Obtained by setting epsilon = 1 
#  The derivative of g
def dg(x):
    c = np.sqrt(8.0)/mu
    return 4.0*x**3 +3.0*c*x**2

#####################################

def h(phi): # solving the integral for e-fold number in phi
    int_manual =  (phi**2-f**2 + mu**4*(phi**-2 - f**-2))/8.0
    return int_manual - N_1                     # N and f are global variables defined below.


def dh(phi):
    x = phi/mu
    return  0.25*mu*(x-x**-3)   

#####################################

for i in range(len(mulist)) :
      mu=mulist[i]
      x0=10.0 
    
      x_end, A, B = newton(g, x0, dg, store = False)
      phi_end= mu*x_end
      phi_end_array.append(phi_end)  
 

phi_list = phi_end_array

for i in range(len(mulist)) :                
    mu=mulist[i]        
    f=phi_list[i]       
    
    for j in range(len(Nlist)) :             

        N_1=Nlist[j]      
        y0 = 0.5*f      
        phi, A, B = newton(h, y0, dh, store = False)
        
        
        phistar = phi/mu
        
        r= f_r(phistar)
        ns= f_ns(phistar)
        f_r_store.append(r)
        f_ns_store.append(ns)

     
       
    plot(f_ns_store,f_r_store, color = 'yellow')     

    f_r_store=[]
    f_ns_store=[]   

r = f_r0(N, n) #phi^2/3
ns = ns_N0(N, n)

r_1 = f_r1(N) #R^2
ns_1 = ns_N1(N)

r_2 = f_r2(N, n_1) #Massive Potential 
ns_2 = ns_N2(N, n_1)

r_3 = f_r3(N, n_2) #phi^4/3
ns_3 = ns_N3(N, n_2)

r_4 = f_r4(N, n_3) #phi^3
ns_4 = ns_N4(N, n_3)

r_5 = f_r5(N, n_4) #phi
ns_5 = ns_N5(N, n_4)

r_6 = f_r6(N, n_5) #Random Potential 
ns_6 = ns_N6(N, n_5)

plot(ns, r, color= 'red', label = '$φ^{2/3}$')
plot(ns_5, r_5, color= 'green', label = 'φ')
plot(ns_6, r_6, color= 'pink', label ='$φ^{1.125}$')
plot(ns_3, r_3, color= 'black', label = '$φ^{4/3}$')
plot(ns_2, r_2, color= 'brown', label = '$φ^2$')
plot(ns_4, r_4, color= 'blue', label = '$φ^3$')
plot(ns_1, r_1, color= 'purple', label = '$R^2$ Inflation')

###################################
alist = np.linspace(0.01, 3500, 500) 
f_r_store_1  = []       
f_ns_store_1 = [] 


def f_r_1(N, a) :
            return 48*a/(4*N**2 + 2*N*sqrt(3*a*(4 + 3*a)) + 3*a)

def f_ns_1(N, a) : 
            return 1 - ((8*N + 6*a + 2*sqrt(3*a*(4 + 3*a)))/(4*N**2 + 2*N*sqrt(3*a*(4 + 3*a)) + 3*a))

for i in range(len(alist)) :                
    a=alist[i]        
          
    
    for j in range(len(Nlist)) :             
        Nlist = np.linspace(50.0, 60.0, 10)
        N=Nlist[j] 



        
        r= f_r_1(N, a)
        ns= f_ns_1(N, a)
        f_r_store_1.append(r)
        f_ns_store_1.append(ns)

    plot(f_ns_store_1, f_r_store_1, color = 'orange')  

    f_r_store_1=[]
    f_ns_store_1=[]
    
plt.axhline(y=0.07, color='brown')
plt.legend()
G.export()
