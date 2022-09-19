import numpy as np
import math
import importPolars as ip
import matplotlib.pyplot as plt
import pandas as pd


def finda(a):
    diff_1 = abs(a[0] - 1 / 3)
    diff_2 = abs(a[1] - 1 / 3)
    diff_3 = abs(a[2] - 1 / 3)

    if diff_1 < diff_2 and diff_1 < diff_3:
        a_opt = a[0]
    elif diff_2 < diff_1 and diff_2 < diff_3:
        a_opt = a[1]
    else:
        a_opt = a[2]

    return a_opt

def calc_u_rel(a, v, da, w, r):
    U_rel = math.sqrt(((1 - a) * v) ** 2 + ((1 + da) * w * r) ** 2)
    return U_rel

def calc_Re(rho, U_rel, L_c, mu):
    Re = (rho * U_rel * L_c) / mu
    return Re

def string_for_save(filename):
    new_filename = filename.replace("/", "_")
    return new_filename

def plot_data(df, x_column, y_column, profile):
    ax = plt.gca()
    df.plot(kind='line', x=x_column, y=y_column, ax=ax, label=profile, ylabel=y_column)
    file_name = string_for_save(f"{y_column}_over_{x_column}_{profile}.png")
    plt.savefig(file_name)
    plt.show()
    ax.clear()

def plot_data_2_graphs(df, x_column, y_column1, y_column2, y_label):
    ax = plt.gca()
    df.plot(kind='line', x=x_column, y=y_column1, ax=ax, label=y_column1, xlabel=x_column, ylabel=y_label)
    df.plot(kind='line', x=x_column, y=y_column2, ax=ax, label=y_column2)
    file_name = string_for_save(f"{y_column1}_and_{y_column2}_over_{x_column}_{profile}.png")
    plt.savefig(file_name)
    ax.clear()



# define the parameters
profile = "S808" #declare the profile name for plotting
R = 0.38  # overall rotorradius (max 0.9m) currently only works uo to 0.6m because of the polar file
r_hub = 0  # hub radius
N = 13  # Number of elements per blade (10...15)
dr = R / N  # distance between element center points

r_vec = np.arange(r_hub + dr / 2, R, dr)  # vector with all the centerpoint location

Z = 2  # blade count

v = 12  # wind speed

lambda_des = 4 * math.pi / Z  - 1.5 # TSR

rho = 1.2  # air density 20°C

mu = 1.47e-5  # dyn. viscosity air 20°C

w = v * lambda_des/R # angular velocity

#initial Re calculation
L_c_init=50e-3 #characteristic length between 40mm and 60mm
a_init =  1/3
da_init= 0
u_rel_init = calc_u_rel(a_init, v, da_init, w, 0.7*R)
Re_init = calc_Re(rho, u_rel_init, L_c_init, mu)  # initial Re guess
tol = 20000 #recommended by institute

# get the polar data for the airfoil
filename = 'S808_Airfoil_11polars.txt'
#TODO: change the file here to a custom file from us
ClFun, CdFun, AoAmax, Relims, AoALims = ip.importPolars(filename,'ashes')
                                                        #'airfoiltools')
                                                        # 'ashes')

df = pd.DataFrame(columns=['r', 'L_c', 'Theta', 'alpha_opt', 'a', 'da', 'Re', 'r/R [m]', 'L_c/R', 'dT', 'dM', 'dM*r']) #create a dataframe to store data

# here starts the loop over all elements
for r in r_vec:
    Re=Re_init
    Re_old=0
    # local tip ratio
    lambda_r = r * lambda_des / R

    # find a and a'
    coeff = [16, -24, 9 - 3 * (lambda_r ** 2), -1 + lambda_r ** 2]
    a_all = np.roots(coeff)
    a = finda(a_all)  # find the a closest to 1/3

    da = (1 - 3 * a) / (4 * a - 1)

    # find flow angle
    Phi = math.atan((1 - a) / ((1 + da) * lambda_r))

    # second loop through all the Re numbers
    while np.absolute(Re - Re_old) >= tol:
        
        Re_old = Re
        # get current Cl, Cd, a_opt
        C_l = ClFun(Re, AoAmax(Re))
        C_d = CdFun(Re, AoAmax(Re))

        # find axial force coefficient
        C_a = C_l * np.cos(Phi) + C_d * np.sin(Phi)

        # find rotational force coefficient
        C_r= C_l * np.sin(Phi) + C_d * np.cos(Phi)

        # find chord length
        L_c = R * (8 * math.pi * a * lambda_r * (np.sin(Phi) ** 2)) / ((1 - a) * Z * C_a * lambda_des)
        # find relative velocity
        U_rel = calc_u_rel(a, v, da, w, r)
        # find new Re
        Re = calc_Re(rho, U_rel, L_c,mu)
        #print('res=', np.absolute(Re - Re_old))
        
    Theta = Phi - AoAmax(Re)
    alpha_opt=AoAmax(Re)
    #Force calculation
    dT = C_a * 0.5 * rho * (U_rel ** 2) * L_c * dr * Z
    dM = C_r * 0.5 * rho * (U_rel ** 2) * L_c * dr * Z
    df.loc[len(df.index)] = [r, L_c, Theta, alpha_opt, a, da, Re, r/R, L_c/R, dT, dM, dM*r] #add values to dataframe
    print('@r/R = ',r/R,', Cord length Lc = ',L_c,', Twist angle Theta = ', Theta,'Re = ',Re)
    
#2.1
print(f"For 2.1:\n"
      f"lambda_des = {lambda_des}\n"
      f"R = {R}m\n"
      f"Z = {Z}\n"
      f"V = {v}m/s\n"
      f"Re_init[at r/R [m] = 0.7] = {Re_init}\n"
      f"Number of elements = {N}")

#TODO: 2.2 is not implemented yet

#2.3
df2 = df[['r/R [m]', 'L_c/R', 'Theta', 'a', 'da', 'alpha_opt']].copy()
df2.to_csv('results_3A_3_S808.csv', index=False, encoding='utf-8')
plot_data(df2,'r/R [m]', 'L_c/R', profile)
plot_data(df2,'r/R [m]', 'Theta', profile)

#2.4
df3 = df[['r/R [m]', 'Re', 'dT', 'dM']].copy()
df3.to_csv('results_3A_4_S808.csv', index=False, encoding='utf-8')
#plot_data_2_graphs(df3,'r/R [m]', 'dT', 'dM', "Force [N]")
M=df["dM*r"].sum()
T=df["dT"].sum()
P=M*w
C_p=P/(0.5*rho*(v**3)*math.pi*(R**2))
C_t=T/(0.5*rho*(v**2)*math.pi*(R**2))
print(f"For 2.4:"
      f"Total power produced: P = {P} W\n"
      f"Power coefficient: C_p = {C_p}\n"
      f"Thrust coefficient: C_t = {C_t}\n")

# Assignment 3B Task 3.1
df4 = df[['r', 'L_c', 'Theta','Re']].copy()
df4.to_csv('results_3B_S808.csv', index=False, encoding='utf-8')    
    




    
    


