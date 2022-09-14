import numpy as np
#import scipy.io
#from scipy import signal
#import matplotlib.pyplot as plt
import math 
import importPolars as ip


def finda(a):
    diff_1 = abs(a[0]-1/3)
    diff_2 = abs(a[1]-1/3)
    diff_3 = abs(a[2]-1/3)
    
    if diff_1 < diff_2 and diff_1 < diff_3:
        a_opt = a[0]
    elif diff_2 < diff_1 and diff_2 < diff_3:
        a_opt = a[1]
    else:
        a_opt = a[2]
    
    return a_opt

# define the parameters, only R and L_init can be varied (lambda_des might aswell but current formula supossedly best)

R = 0.9             # overall rotorradius (max 0.9m)
r_hub = 0           # hub radius
N = 10              # Number of elements per blade (10...15)
dr = R/N            # distance between element center points

r_vec = np.arange(r_hub+dr/2,R,dr) # vector with all the centerpoint location

Z = 2               # blade count

v = 12               # wind speed

lambda_des = 4*math.pi/Z      # TSR

rho = 1.2         # air density 20°C

mu = 1.47e-5        # dyn. viscosity air 20°C

w = lambda_des*v/R     # rotational speed

W_init = math.sqrt(((1-1/3)*v)**2 + ((1+0)*w*0.7*R)**2)
L_init = 0.05
Re_init = rho*W_init*L_init/mu  # initial Re guess
Re = Re_init

Re_old = 0              # set for iteration loop criterion
tol = 20000    

# get the polar data for the airfoil
filename = 'S826_8polars18.txt'
ClFun, CdFun, AoAmax, Relims, AoALims = ip.importPolars(filename,'ashes')     

# here starts the loop over all elements
for r in r_vec:
    
    # local tip ratio
    lambda_r = r*lambda_des/R
    
    
    # find a and a'
    coeff = [16,-24,9-3*(lambda_r**2),-1+lambda_r**2]
    a_all = np.roots(coeff)
    a = finda(a_all)            # find the a closest to 1/3
    
    da = (1-3*a)/(4*a-1)
    
    
    # find flow angle   
    Phi = math.atan((1-a)/((1+da)*lambda_r))

    # second loop through all the Re numbers 
    while np.absolute(Re - Re_old) >= tol:
        
        Re_old = Re
        
        # get current Cl, Cd, a_opt
        C_l = ClFun(Re,AoAmax(Re))
        C_d = CdFun(Re,AoAmax(Re))
        
        # find thrust coefficient
        C_a = C_l*np.cos(Phi) + C_d*np.sin(Phi)
        
        # find chord length
        L_c = R*(8*math.pi*a*lambda_r*( np.sin(Phi)**2 )) / ((1-a)*Z*C_a*lambda_des)
        
        # find relative velocity
        U_rel = math.sqrt( ((1-a)*v)**2 + ((1+da)*w*r)**2 )
        
        # find new Re
        Re = (rho*U_rel*L_c) / mu
        
    Theta = Phi - AoAmax(Re)
    
    print('@r = ',r,', Cord length Lc = ',L_c,', Twist angle Theta = ', Theta)


