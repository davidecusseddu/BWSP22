# This script is attached to 'Numerical investigations of the bulk-surface wave pinning model' by Cusseddu and Madzvamuse 2022
# It solves the BSWP model presented in that work
# For more details see also Section 5 in Cusseddu et al., JTB 19 (doi = 10.1016/j.jtbi.2018.09.008)
# Main file: BSWPmodel.py 

import numpy as np
import os

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

# TEMPORAL DISCRETISATION
t0 = 0.0 # initial time
Tf = 15005 # final time
tau = 1e-2 # discretisation step

# SPATIAL DISCRETISATION - INPUT MESH FROM KEYBOARD
print('Select domain and mesh:\n')
print('1. Oblate spheroid')
print('2. Ellipse')
print('3. Egg')
print('4. Non-convex domain\n')
selected_mesh = input('type the number of the mesh:\n') 

# DEFINITION OF THE INITIAL CONDITION FOR a
print('The code will automatically set the initial conditions as in the paper. Access the code to change them\n')

if selected_mesh == '1':

    str_mesh = 'oblate'
    
    a_p0 = 1
    x0 = 0
    y0 = 0.85
    z0 = 0
    invsigma2_x0 = 10
    invsigma2_y0 = 10
    invsigma2_z0 = 10
    m0 = 1 # m0 = M0/|Omega|

    ic_a = 'a_p0*(exp(-invsigma2_x0*pow(x[0]-x0,2)-invsigma2_y0*pow(x[1]-y0,2)-invsigma2_z0*pow(x[2]-z0,2) ))'
    
elif selected_mesh == '2':

    str_mesh = 'ellipse'

    ic_a_top = 'a_p0*(exp(-invsigma2_x0*pow(x[0]-x0,2)-invsigma2_y0*pow(x[1]-y0,2) )) '
    ic_a_bottom = 'a_p1*(exp(-invsigma2_x1*pow(x[0]-x1,2)-invsigma2_y1*pow(x[1]-y1,2)))'
    ic_a = ic_a_top + ' + ' + ic_a_bottom

    a_p0 = 1
    x0 = 0
    y0 = 0.85
    z0 = 0
    invsigma2_x0 = 10
    invsigma2_y0 = 10
    
    a_p1 = 1
    x1 = 0
    y1 = -0.85
    z1 = 0
    invsigma2_x1 = 10
    invsigma2_y1 = 10

    m0 = 0.85 # m0 = M0/|Omega|

    
elif selected_mesh == '3':

    refinement = input('Please specify the level of refinment: 0 (less refined), 1, 2 (more refined) \n = ')
    if refinement == '1':
        str_mesh = 'egg'
    else:
        str_mesh = 'egg_refined'

    ic_a_top = 'a_p0*(exp(-invsigma2_x0*pow(x[0]-x0,2)-invsigma2_y0*pow(x[1]-y0,2)-invsigma2_z0*pow(x[2]-z0,2)))'
    ic_a_bottom = 'a_p1*(exp(-invsigma2_x1*pow(x[0]-x1,2)-invsigma2_y1*pow(x[1]-y1,2)-invsigma2_z1*pow(x[2]-z1,2)))'
    ic_a = ic_a_top + ' + ' + ic_a_bottom
    a_p0 = 1
    x0 = 0
    y0 = 1.05
    z0 = 0
    invsigma2_x0 = 10
    invsigma2_y0 = 10
    invsigma2_z0 = 10
    a_p1 = 1
    x1 = 0
    y1 = -0.95
    z1 = 0
    invsigma2_x1 = 10
    invsigma2_y1 = 10
    invsigma2_z1 = 10

    m0 = 1 # m0 = M0/|Omega|
    
elif selected_mesh == '4':

    str_mesh = 'cones'
    ic_a =  'a_p0*(exp(-invsigma2_x0*pow(x[0]-x0,2)-invsigma2_y0*pow(x[1]-y0,2)-invsigma2_z0*pow(x[2]-z0,2) ))'
    
    a_p0 = 1
    x0 = 0
    y0 = 2.6
    z0 = 0
    invsigma2_x0 = 0
    invsigma2_y0 = 15
    invsigma2_z0 = 0
    
    m0 = 1.2 # m0 = M0/|Omega|

else:
    raise ValueError("No mesh has been selected")

# SELECT BETWEEN BSWP MODEL OR SURFACE MODEL (REDUCED MODEL) - INPUT FROM KEYBOARD
model = ''
print('')
while model != 'BS' and model != 'S':
    model = input('If you want to solve the BSWP model digit BS, if the reduced model digit S:\n model = ')
    if model=='bs':
        model = 'BS'
    elif model == 's' or model =='S':
        model = 'S'
        Db = '+inf'
print('\n')

# Set the value of the bulk diffusion coefficient
if model == 'BS':
    Db = float(input('Please specify a bulk diffusion coefficient\n Db = '))    



    
    
# CREATE STRING FOR EXPORTING IC PARAMETERS TO TEXT FILE

initial_condition_string = '\n ic_a = ' + ic_a
initial_condition_string = initial_condition_string + '\n a_p0 = ' + str(a_p0)
initial_condition_string = initial_condition_string + '\n x0 = ' + str(x0)
initial_condition_string = initial_condition_string + '\n y0 = ' + str(y0)

if selected_mesh != '2':
    initial_condition_string = initial_condition_string + '\n z0 = ' + str(z0)
    
initial_condition_string = initial_condition_string + '\n invsigma2_x0 = ' + str(invsigma2_x0)
initial_condition_string = initial_condition_string + '\n invsigma2_y0 = ' + str(invsigma2_y0)

if selected_mesh != '2':
    initial_condition_string = initial_condition_string + '\n invsigma2_z0 = ' + str(invsigma2_z0)

if selected_mesh == '2' or selected_mesh == '3':
    initial_condition_string = initial_condition_string + '\n a_p1 = ' + str(a_p1)
    initial_condition_string = initial_condition_string + '\n x1 = ' + str(x1)
    initial_condition_string = initial_condition_string + '\n y1 = ' + str(y1)

    if selected_mesh != '2':
        initial_condition_string = initial_condition_string + '\n z1 = ' + str(z1)

    initial_condition_string = initial_condition_string + '\n invsigma2_x1 = ' + str(invsigma2_x1)
    initial_condition_string = initial_condition_string + '\n invsigma2_y1 = ' + str(invsigma2_y1)

    if selected_mesh != '2':
        initial_condition_string = initial_condition_string + '\n invsigma2_z1 = ' + str(invsigma2_z1)
    
initial_condition_string = initial_condition_string + '\n m0 = ' + str(m0)
    
print(initial_condition_string)



# CREATE STRING FOR STORING SOLUTIONS AND DATA
filename = 'mesh-' + str_mesh + '-' + model + '-k0-' + str(k0) + '-gamma-'+str(gamma) + '-K-'+str(K) + '-beta-' + str(beta) + '-m0-' + str(m0) + '-Da-'+str(Da) + '-Db-'+str(Db)

# DEFINITION OF SOLUTION SAVING TIME POINTS  
savingtimes = [0.0, 1.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 450.0, 550.0, 650.0, 800.0, 1000.0, 1200.0, 1500.0, 2000.0, 2500.0, 3000.0, 4000.0, 5000.0, 8000.0, 11000.0, 15000]
    
print('--------------\n')
print('Resolution of the ', model, ' model\n')
print('The solution will be saved at the timepoints: ')
print([int(x) for x in savingtimes])
print('change this choice in parameters.py\n')
print('The name associated to the solution is')
print(filename)
print('--------------\n')

##############################################
### Create folder for storing system solutions and data:
from datetime import datetime
# datetime object containing current date and time
now = datetime.now()
print("now =", now)

folder = 'output' + str(now) + '/'
os.makedirs(folder[:-1])#,exist_ok=True)
 

parameters_and_data = 'Mesh = ' + str_mesh + \
'\n model = ' + model + \
'\n\n PARAMETERS ' + \
'\n\n Reaction parameters: ' + \
'\n k0 = ' + str(k0) + \
'\n gamma = ' + str(gamma) +\
'\n K = ' + str(K) +\
'\n beta = ' + str(beta) +\
'\n\n Diffusion parameters: ' + \
'\n Da = ' + str(Da) +\
'\n Db = ' + str(Db) +\
'\n\n Time parameters: ' + \
'\n t0 = ' + str(t0) +\
'\n Tf = ' + str(Tf) +\
'\n tau = ' + str(tau) +\
'\n\n Initial conditions: \n' + initial_condition_string
    
    
print('---')
print('Data are exported in ')
with open(folder + filename + '-parameters_and_data.txt', "w") as text_file:
    text_file.write(parameters_and_data)
print(folder + '/' + filename + '-parameters_and_data.txt')
print('---')



