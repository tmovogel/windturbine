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

def check_L(theta, Lc):
    Lmax = 0.098
    Lmin = 0.04
    L = Lc*np.cos(np.radians(theta))
    if L > Lmax:
        Lc = abs(Lmax/np.cos(np.radians(theta)))
    elif L < Lmin:
        Lc = Lmin
    return Lc



# define the parameters
profile = "S809" #declare the profile name for plotting
filename = 'airfoils/' + profile + '_polars.txt' # get the polar data for the airfoil
R = 0.45  # overall rotorradius (max 0.9m) currently only works uo to 0.6m because of the polar file
r_hub =  0.03 # hub radius
N = 15  # Number of elements per blade (10...15)
#dr = R / (N)  # distance between element center points

#r_vec = np.arange(r_hub + dr / 2, R, dr)  # vector with all the centerpoint location
r_vec = np.linspace(r_hub, R, num = N + 1)
dr = r_vec[2]-r_vec[1];
r_vec = r_vec[:N]+dr/2

Z = 2  # blade count

v = 12  # wind speed

lambda_des = 4 * math.pi / Z  -0.2 # TSR

rho = 1.2  # air density 20°C

mu = 1.47e-5  # dyn. viscosity air 20°C



w = v * lambda_des/R # angular velocity

#initial Re calculation
L_c_init = 5e-3 #characteristic length between 40mm and 60mm
a_init =  1/3
da_init= 0
u_rel_init = calc_u_rel(a_init, v, da_init, w, 0.7*R)
Re_init = calc_Re(rho, u_rel_init, L_c_init, mu)  # initial Re guess
tol = 15000 #recommended by institute


#TODO: change the file here to a custom file from us
ClFun, CdFun, AoAmax, Relims, AoALims = ip.importPolars(filename,'ashes')
                                                        #'airfoiltools')
                                                        # 'ashes')

df = pd.DataFrame(columns=['r', 'Theta', 'L_c', 'alpha_opt', 'a', 'da', 'Re', 'r/R [m]', 'L_c/R', 'dT', 'dM', 'dM*r', 'r-r_hub', 'dL']) #create a dataframe to store data

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
    
    it = 0
    # second loop through all the Re numbers
    while  np.absolute(Re - Re_old) > tol:
        it = it + 1
        #print(it)
        
        Re_old = Re
        # get current Cl, Cd, a_opt
        C_l = ClFun(Re, AoAmax(Re))
        C_d = CdFun(Re, AoAmax(Re))

        # find axial force coefficient
        C_a = C_l * np.cos(Phi) + C_d * np.sin(Phi)

        # find rotational force coefficient
        C_r= C_l * np.sin(Phi) + C_d * np.cos(Phi)

        # find chord length
        Lc = R * (8 * math.pi * a * lambda_r * (np.sin(Phi) ** 2)) / ((1 - a) * Z * C_a * lambda_des)
        # find relative velocity
        U_rel = calc_u_rel(a, v, da, w, r)
        # find new Re
        Re = abs(calc_Re(rho, U_rel, Lc,mu))
        #print('res=', np.absolute(Re - Re_old))
        
        
        
    Theta = np.degrees(Phi) - AoAmax(Re)
    alpha_opt = AoAmax(Re)
    L_c = check_L(Theta, Lc)
    dL = abs(L_c-Lc)
    
    #Force calculation
    dT = C_a * 0.5 * rho * (U_rel ** 2) * L_c * dr * Z
    dM = C_r * 0.5 * rho * (U_rel ** 2) * L_c * dr * Z
    df.loc[len(df.index)] = [r, L_c, Theta, alpha_opt, a, da, Re, r/R, L_c/R, dT, dM, dM*r, r-r_hub, dL] #add values to dataframe
    print('Lc = ',L_c,'/',Lc, '   dL = ', dL)
    
#2.1
print(f""
      f"For 2.1:\n"
      f"lambda_des = {lambda_des}\n"
      f"R = {R}m\n"
      f"Z = {Z}\n"
      f"V = {v}m/s\n"
      f"Re_init[at r/R [m] = 0.7] = {Re_init}\n"
      f"Number of elements = {N}")


# #2.3
# df2 = df[['r/R [m]', 'L_c/R', 'Theta', 'a', 'da', 'alpha_opt']].copy()
# df2.to_csv('results_3A_3_S808.csv', index=False, encoding='utf-8')
# # plot_data(df2,'r/R [m]', 'L_c/R', profile)
# # plot_data(df2,'r/R [m]', 'Theta', profile)
#
#2.4
df3 = df[['r/R [m]', 'Re', 'dT', 'dM']].copy()
#df3.to_csv('results_3A_4_S807.csv', index=False, encoding='utf-8')
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
df4 = df[['r','Theta', 'L_c' ,'Re', 'r-r_hub', 'alpha_opt']].copy()
#df["L_c"] = np.where(df["L_c"] < 0.04, 0.04, df['L_c'])
file = 'results_' + profile + '.csv'
df4.to_csv(file, index=False, encoding='utf-8')
