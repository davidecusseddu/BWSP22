# This script is attached to 'Numerical investigations on the bulk-surface wave pinning model' by Cusseddu and Madzvamuse 2022
# It solves the BSWP model presented in that work
# For more details see also Section 5 in Cusseddu et al., JTB 19 (doi = 10.1016/j.jtbi.2018.09.008)
# Main file: BSWPmodel.py 

import numpy as np

# DEFINITION OF THE SYSTEM PARAMETERS / DISCRETISATION / STORING DATA

# REACTION FUNCTION
k0 = 0.1 # basal activation rate
gamma = 1.0 # positive feedback activation rate
K = 0.3 # concentration parameter positive feedback
beta = 1.0 # basal inactivation rate

def f1(a): # activation term to multiply by b
    return (k0+gamma*(a**2)/(K**2+a**2))
def f(a,b): # full reaction function
    return b*f1(a)-beta*a

# DIFFUSION
Da = 0.005 # surface diffusion parameter
Db = 0.5 # bulk diffusion parameter

# TEMPORAL DISCRETISATION
t0 = 0.0 # initial time
Tf = 15005 # final time
tau = 1e-2 # discretisation step

# SPATIAL DISCRETISATION - INPUT MESH FROM KEYBOARD
print('Select domain and mesh:\n')
print('1. For the oblate spheroid type oblate3d')
print('2. For its 2D version type oblate2d')
print('3. For the non-convex domain type cones\n')
str_mesh = input('type the name of the mesh:\n') 

# SELECT BETWEEN BSWP MODEL OR SURFACE MODEL (REDUCED MODEL) - INPUT FROM KEYBOARD
model = ''
print('')
while model != 'BS' and model != 'S':
    model = input('If you want to solve the BSWP model digit BS, if the reduced model digit S:\n model = ')
    if model=='bs':
        model = 'BS'
    elif model == 's':
        model = 'S'
        Db = '+inf'
print('\n')

# INITIAL CONDITION FOR a
if str_mesh == 'oblate2d':
    ic_a = 'C_up*(exp(-invsigma2_x*pow(x[0]-x0,2)-invsigma2_y*pow(x[1]-y0,2)))'
else:
    ic_a = 'C_up*(exp(-invsigma2_x*pow(x[0]-x0,2)-invsigma2_y*pow(x[1]-y0,2)-invsigma2_z*pow(x[2]-z0,2)))'
C_up = 1 # height gaussian peak
invsigma2_x=10 # bell width parameter x
invsigma2_y=10 # bell width parameter y
invsigma2_z=10 # bell width parameter z
x0 = 0 # center of gaussian peak
y0 = 0.85 # center of gaussian peak
z0 = 0 # center of gaussian peak


# CONSERVATION OF MASS
m0 = 1. # m0 = M0/|Omega|


# CREATE STRING FOR STORING SOLUTIONS AND DATA
filename = 'mesh-' + str_mesh + '-' + model + '-k0-' + str(k0) + '-gamma-'+str(gamma) + '-K-'+str(K) + '-beta-' + str(beta) + '-m0-' + str(m0) + '-Da-'+str(Da) + '-Db-'+str(Db)

# DEFINITION OF SOLUTION SAVING TIME POINTS  
savingtimes = [0.0, 1.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 450.0, 550.0, 650.0, 800.0, 1000.0, 1200.0, 1500.0, 2000.0, 2500.0, 3000.0, 4000.0, 5000.0, 8000.0, 11000.0, 15000]

# METHOD FOR SOLVING THE LINEAR SYSTEMS 
solver_method = 'cg' # CONJUGATE GRADIENT METHOD
preconditioner = 'sor' 

    
print('--------------\n')
print('Resolution of the ', model, ' model\n')
print('The solution will be saved at the timepoints: ')
print([int(x) for x in savingtimes])
print('change this choice in parameters.py\n')
print('The name associated to the solution is')
print(filename)
print('--------------\n')

