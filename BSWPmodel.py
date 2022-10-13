# This script is attached to 'Numerical investigations of the bulk-surface wave pinning model' by Cusseddu and Madzvamuse 2022
# It solves the BSWP model presented in that work
# For more details see also Section 5 in Cusseddu et al., JTB 19 (doi = 10.1016/j.jtbi.2018.09.008)
# Matrices are defined using same notation as Cusseddu et al. JTB 19

#####################################
# IMPORT LIBRARIES:
from fenics import *
from dolfin import *

from math import *
import numpy as np


##############################################
### Import parameter values from parameter_values.py:
from parameter_values import *
parameters["form_compiler"]["quadrature_degree"] = 5


##############################################
### Import the mesh and extract boundary mesh:
mesh = 'mesh/'+str_mesh + '.xml' # The mesh on Omega is given in input from keyboard (see parameter_values.py)
domain_mesh = Mesh(mesh)

if selected_mesh == '3' and refinement == '2':
    domain_mesh = refine(domain_mesh)

dX = dx(domain_mesh)
boundary_mesh = BoundaryMesh(domain_mesh, 'exterior') # The boundary mesh is extracted from bulk mesh
dx = dx(boundary_mesh)



##############################################
### Definition of the function spaces over each mesh:
if model == 'BS':
    V_Bulk = FunctionSpace(domain_mesh, 'P', 1) # not needed for the reduced model    
V_Boundary = FunctionSpace(boundary_mesh,'P', 1)



##############################################
### Extraction of the relationship between the vertices of the two meshes:
if model == 'BS': 
    # Mapping vertex - dof between domain_mesh and boundary_mesh:
    dofb_vboundary = np.array(dof_to_vertex_map(V_Boundary), dtype=int)
    vboundary_v = boundary_mesh.entity_map(0).array()
    vBulk_dof = np.array(vertex_to_dof_map(V_Bulk), dtype=int)
    parameters['allow_extrapolation'] = True

    

##############################################
### Extract mesh properties and save in a text file:
bulk_volume = assemble(Constant(1.0)*dX) # Calculate volume of Omega |Omega|
surface_area = assemble(Constant(1.0)*dx) # Calculates Gamma surface area |Gamma|
coordinate_Omega=domain_mesh.coordinates()
coordinate_Boundary=boundary_mesh.coordinates()
cellsOmega=domain_mesh.cells()
cellspartialOmega=boundary_mesh.cells()
mesh_str = 'Mesh = ' + str_mesh + \
'\n n vertices in bulk mesh = ' + str(domain_mesh.num_vertices()) + \
'\n n vertices on boundary mesh = ' + str(boundary_mesh.num_vertices()) + \
'\n ( n vertices in bulk mesh, dimension ) = ' + str(np.shape(coordinate_Omega)) + \
'\n ( n vertices in boundary mesh, dimension ) = ' + str(np.shape(coordinate_Boundary)) + \
'\n ( cells in bulk mesh, number of vertices per cell ) = ' + str(np.shape(cellsOmega)) + \
'\n ( cells in boundary mesh, number of vertices per cell ) =' + str(np.shape(cellspartialOmega)) + \
'\n Bulk mesh: (hmax, hmin) = ' + str([domain_mesh.hmax(), domain_mesh.hmin()]) + \
'\n hmax-hmin = ' + str(domain_mesh.hmax()-domain_mesh.hmin()) + \
'\n Boundary mesh: (hmax, hmin) = ' + str([boundary_mesh.hmax(), boundary_mesh.hmin()]) + \
'\n hmax-hmin = ' + str(boundary_mesh.hmax()-boundary_mesh.hmin()) + \
'\n bulk_volume = ' + str(bulk_volume) + \
'\n surface_area = ' + str(surface_area)
print('---')
print('Mesh properties are exported in ')
with open(folder + filename + '-mesh_data.txt', "w") as text_file:
    text_file.write(mesh_str)
print(folder + '/' + filename + '-mesh_data.txt')
print('---')



##############################################
### Setting the initial conditions:
# Initial condition for a
print('Setting up the initial conditions...\n')

if selected_mesh == '1' or selected_mesh =='4':
    a_0 = Expression(ic_a, degree = 3, a_p0=a_p0, invsigma2_x0=invsigma2_x0, x0=x0, invsigma2_y0=invsigma2_y0, y0=y0, invsigma2_z0=invsigma2_z0, z0=z0)    

elif selected_mesh == '2':
    a_T = Expression(ic_a_top, degree = 3, a_p0=a_p0, invsigma2_x0=invsigma2_x0, x0=x0, invsigma2_y0=invsigma2_y0, y0=y0)
    a_B = Expression(ic_a_bottom, degree = 3, a_p1=a_p1, invsigma2_x1=invsigma2_x1, x1=x1, invsigma2_y1=invsigma2_y1, y1=y1)
    a_T_n = interpolate(a_T, V_Boundary)
    a_B_n = interpolate(a_B, V_Boundary)
    print('mass of a top = ', assemble(a_T_n*dx))
    print('mass of a bottom = ', assemble(a_B_n*dx))
    a_0 = Expression(ic_a, degree = 3, a_p0=a_p0, invsigma2_x0=invsigma2_x0, x0=x0, invsigma2_y0=invsigma2_y0, y0=y0, a_p1=a_p1, invsigma2_x1=invsigma2_x1, x1=x1, invsigma2_y1=invsigma2_y1, y1=y1)

elif selected_mesh == '3':
    a_T = Expression(ic_a_top, degree = 3, a_p0=a_p0, invsigma2_x0=invsigma2_x0, x0=x0, invsigma2_y0=invsigma2_y0, y0=y0, invsigma2_z0=invsigma2_z0, z0=z0)
    a_B = Expression(ic_a_bottom, degree = 3, a_p1=a_p1, invsigma2_x1=invsigma2_x1, x1=x1, invsigma2_y1=invsigma2_y1, y1=y1, invsigma2_z1=invsigma2_z1, z1=z1)
    a_T_n = interpolate(a_T, V_Boundary)
    a_B_n = interpolate(a_B, V_Boundary)
    print('mass of a top = ', assemble(a_T_n*dx))
    print('mass of a bottom = ', assemble(a_B_n*dx))
    a_0 = Expression(ic_a, degree = 3, a_p0=a_p0, invsigma2_x0=invsigma2_x0, x0=x0, invsigma2_y0=invsigma2_y0, y0=y0, invsigma2_z0=invsigma2_z0, z0=z0, a_p1=a_p1, invsigma2_x1=invsigma2_x1, x1=x1, invsigma2_y1=invsigma2_y1, y1=y1, invsigma2_z1=invsigma2_z1, z1=z1)


a_n = interpolate(a_0, V_Boundary)

# Initial condition for b
b_homogeneous_n = m0 - assemble(a_n*dx)/bulk_volume
if model == 'BS':
    b_n = interpolate(Constant(b_homogeneous_n), V_Bulk)

    b_boundary = interpolate(b_n, V_Boundary)     # Restriction of b_n to the boundary
    
    # For compatibility issues, we extend a over the bulk by creating the following
    a_bulk_n = Function(V_Bulk)
    array = a_bulk_n.vector().get_local() # Extraction of the ordered values of a_bulk_n in array
    array[vBulk_dof[vboundary_v[dofb_vboundary]]] = a_n.vector().get_local() # only the values on the boundary are replaced by the values of a_n
    a_bulk_n.vector()[:] = array # now a_bulk_n corresponds to a_n over Gamma and it is zero elsewhere
    
    a_bulk_p = Function(V_Bulk)  # Extension of the predicted a_p to the function space V_Bulk



##############################################
### Definition of the trial and test functions:
if model == 'BS':
    b = TrialFunction(V_Bulk)
    v = TestFunction(V_Bulk)
a = TrialFunction(V_Boundary)
w = TestFunction(V_Boundary)



##############################################
### Assemblying the systems:

# Surface system - LHS
M_Gamma = assemble(a*w*dx) #mass matrix
K_Gamma = assemble(inner(grad(a), grad(w))*dx) #stiffness matrix
A_Gamma_P = M_Gamma + tau*Da*K_Gamma # Matrix of the coefficients predictor step (P)
A_Gamma_C = M_Gamma + 0.5*tau*Da*K_Gamma # Matrix of the coefficients corrector step (C)

# Definition of the surface unknowns
a_p = Function(V_Boundary) # predictor for a
a = Function(V_Boundary) # corrector for a

if model == 'BS':    
    # Bulk system - LHS:
    M_Omega = assemble(b*v*dX) #mass matrix
    K_Omega = assemble(inner(grad(b), grad(v))*dX) #stiffness matrix
    A_Omega_fixed = M_Omega + 0.5*tau*Db*K_Omega # Stationary part of the matrix of the coefficients
    RHS_Omega_fixed = M_Omega - 0.5*tau*Db*K_Omega # Stationary part in the RHS (it needs to be multiplied by b_n at each time)
    
    # Definition of the bulk unknown
    b = Function(V_Bulk)



##############################################
# Creation of pvd files for storing the solutions:
Solution_a = File(folder + '/' + filename + '-Solution_a.pvd')
if model == 'BS':
    Solution_b = File(folder + '/' + filename + '-Solution_b.pvd')    



##############################################
# Vectors a_p, b and a of the solution values over the dofs:
A_p = a_p.vector()
if model == 'BS':
	B = b.vector()
A = a.vector()




#####################################
# Inizialization time, time steps, vector for ||a - a_n||_2
t = t0
time=[] # List where time steps will be stored
l2norm_Delta_a = [] # List where the L2-norm of a - a_n will be stored at each time
mass_a = []


i=0
index_saving = 0 #13 ## 0
lensavingtimes = len(savingtimes)

m_a = assemble(a_n*dx)
mass_a.append(m_a)
if model == 'BS':
    mass_tot = [m_a+ assemble(b_n*dX)] # List where the integral of b will be stored at each time


#####################################
### SOLVING THE SYSTEM AT EACH TIME t
print('...Solving the system...')
timer_solution = Timer('eval') # Starting the counting of the computational time

while t<Tf:
    t += tau # Increases the time step

    #####################################
    ### STORING THE SOLUTION
    if index_saving < lensavingtimes and t >= savingtimes[index_saving]:
        Solution_a << (a_n, t-tau)
        np.savetxt(folder + '/' + filename + '-last_solution_a_t='+ str(t-tau) + '.txt', a_n.vector())
        np.savetxt(folder + '/' + filename + '-mass_a-t=' + str(t-tau) + '.txt', mass_a)
        mass_a = []
        if model == 'BS':
            Solution_b << (b_n, t-tau)
            np.savetxt(folder + '/' + filename + '-last_solution_b_t='+ str(t-tau) + '.txt', b_n.vector())
            np.savetxt(folder + '/' + filename + '-mass_tot-t=' + str(t-tau) + '.txt', mass_tot)
            mass_tot = []
        index_saving += 1
        print('Solution saved at time t = ', int(t),'. Time required ', timer_solution.stop(), ' s') 
        timer_solution = Timer('eval') 
        timer_solution.start() # Restarting the counting of the computational time



    #####################################
    ### SOLVING THE FIRST SYSTEM - predictor for a EQ (41) of Cusseddu et al., JTB 2019
    # Assemblying the right hand side
    M_Gamma_a_n = M_Gamma*a_n.vector() 
    if model == 'BS':
        F_n = assemble(f(a_n,b_boundary)*w*dx)
        RHS_Gamma_P =  tau*F_n + M_Gamma_a_n 
    else:
        F1_n = assemble(f1(a_n)*w*dx)
        RHS_Gamma_P =  tau*b_homogeneous_n*F1_n + (1-tau*beta)*M_Gamma_a_n

    # Solving the first system using solver_method with preconditioner
    solve(A_Gamma_P, A_p, RHS_Gamma_P) ## , solver_method, preconditioner) # it returns a_p, predictor for a



    #####################################
    ### SOLVING THE SECOND SYSTEM - Solution b EQ (42) of of Cusseddu et al., JTB 2019
    # Assemblyng the right hand side
    if model == 'BS':

        # Extension of the predicted a to the bulk function space (as described above)
        array = a_bulk_p.vector().get_local()
        array[vBulk_dof[vboundary_v[dofb_vboundary]]] = a_p.vector().get_local()
        a_bulk_p.vector()[:] = array
        
        # Assemblying the left hand side
        bulkfunction = TrialFunction(V_Bulk) # I define this function for creating the matrix G_p
        G_p = assemble(f1(a_bulk_p)*bulkfunction*v*ds)
        A_Omega = A_Omega_fixed + 0.5*tau*G_p # The coefficient matrix is now completed

        # Assemblying the right hand side
        G_n = assemble(f1(a_bulk_n)*b_n*v*ds)
        H_a_p = assemble(a_bulk_p*v*ds)
        H_a_n = assemble(a_bulk_n*v*ds)
        RHS_Omega = RHS_Omega_fixed*b_n.vector() - 0.5*tau*G_n + 0.5*tau*beta*(H_a_p+H_a_n)

        # Solving the second system using solver_method with preconditioner
        solve(A_Omega, B, RHS_Omega) ## , solver_method, preconditioner) # it returns b
        
    	# For compatibility issues, the function b is restricted to V_Boundary
        b_boundary_new = interpolate(b,V_Boundary)

    else:

        # Calculate the value of b using conservation of total mass
        b_homogeneous_p = m0 - assemble(a_p*dx)/bulk_volume




    #####################################
    ### SOLVING THE THIRD SYSTEM - Solution b EQ (43) of of Cusseddu et al., JTB 2019
    # Assemblyng the right hand side
    Ms_p = M_Gamma*a_p.vector()
    if model == 'BS':
        F_p = assemble(f(a_p,b_boundary_new)*w*dx)
        RHS_Gamma_C = M_Gamma_a_n - 0.5*tau*Da*K_Gamma*a_n.vector() + 0.5*tau*(F_n + F_p)
    else:
        F1_p = assemble(f1(a_p)*w*dx)
        RHS_Gamma_C = (1-0.5*tau*beta)*M_Gamma_a_n - 0.5*tau*beta*Ms_p - 0.5*tau*Da*K_Gamma*a_n.vector() + 0.5*tau*(b_homogeneous_n*F1_n + b_homogeneous_p*F1_p)
    
    # Solving the third system using solver_method with preconditioner
    solve(A_Gamma_C, A, RHS_Gamma_C) ## solver_method, preconditioner) # it returns a, corrector for a




    #####################################
    ### L2-norm of a-a_n is calculated and added to its list. Also time is saved, and the integral of a
    time.append(t)
    l2norm_Delta_a.append(sqrt(assemble((a-a_n)**2*dx))/tau)
    m_a = assemble(a*dx)
    mass_a.append(m_a)


    #####################################
    ### Update solutions:
    a_n.assign(a)
    if model == 'BS':
        b_n.assign(b)
        b_boundary.assign(b_boundary_new)
    
        a_bulk_n = Function(V_Bulk)
        array = a_bulk_n.vector().get_local()
        array[vBulk_dof[vboundary_v[dofb_vboundary]]] = a_n.vector().get_local()
        a_bulk_n.vector()[:] = array

        mass_tot.append(m_a + assemble(b*dX))
        
    else:
        b_homogeneous_n = m0 - assemble(a_n*dx)/bulk_volume



#####################################
# Saving the l2norm of a -a_n, as well as the last calculated solutions, in text filea:
np.savetxt(folder + '/' + filename + '-l2norm_Delta_a.txt', l2norm_Delta_a)
np.savetxt(folder + '/' + filename + '-last_solution_a_t='+ str(t) + '.txt', a_n.vector())
np.savetxt(folder + '/' + filename + '-mass_a-t=' + str(t) + '.txt', mass_a)

if model == 'BS':
    np.savetxt(folder + '/' + filename + '-last_solution_b_t='+ str(t) + '.txt', b_n.vector())
    np.savetxt(folder + '/' + filename + '-mass_tot-t=' + str(t) + '.txt', mass_tot)




print('....')
print('Process finished.')
print('Results saved in the folder ', folder)
print('reference string:', filename)
