# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.optimize import minimize

# def flory_huggins_deriv(chi, phi1):
#     return  np.log(phi1) - np.log(1 - phi1) +  chi*(1 - 2*phi1)

# def objective_function(x, chi):
#     phi1A = x[0]
#     phi1B = x[1]
#     muA = flory_huggins_deriv(chi, phi1A)
#     muB = flory_huggins_deriv(chi, phi1B)
#     return (muA - muB)**2

# def composition_solver(chi):
#     bounds = [(0,1),(0,1)]
#     result = minimize(lambda x:objective_function(x,chi),[1e-3,1-1e-3],bounds=bounds,method="trust-ncg")
#     print(result.fun)
#     return result.x

# vals = composition_solver(4.5)
# print(vals)

# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.optimize import minimize

# def flory_huggins_deriv(chi, phi1):
#     return np.log(phi1) - np.log(1 - phi1) + chi * (1 - 2 * phi1)

# def objective_function(x, chi):
#     phi1A = x[0]
#     phi1B = x[1]
#     muA = flory_huggins_deriv(chi, phi1A)
#     muB = flory_huggins_deriv(chi, phi1B)
#     return (muA - muB)**2

# def jacobian(x, chi):
#     phi1A = x[0]
#     phi1B = x[1]
    
#     dmuA_dphi1A = (1 / phi1A) + chi * (-2) - (1 / (1 - phi1A))
#     dmuB_dphi1B = (1 / phi1B) + chi * (-2) - (1 / (1 - phi1B))
    
#     dobj_dphi1A = 2 * (flory_huggins_deriv(chi, phi1A) - flory_huggins_deriv(chi, phi1B)) * dmuA_dphi1A
#     dobj_dphi1B = 2 * (flory_huggins_deriv(chi, phi1A) - flory_huggins_deriv(chi, phi1B)) * (-dmuB_dphi1B)
    
#     return np.array([dobj_dphi1A, dobj_dphi1B])

# def composition_solver(chi):
#     bounds = [(1e-12, 1 - 1e-12), (1e-12, 1 - 1e-12)]
#     result = minimize(lambda x: objective_function(x, chi), [1e-3, 1 - 1e-3], 
#                       bounds=bounds, method="trust-constr", jac=lambda x: jacobian(x, chi))
#     print(result.fun)
#     return result.x

# vals = composition_solver(2.5)
# print(vals)


# def solve_equilibrium_compositions(Nchi, N1_N2):
#     initial_guesses = [0.01, 0.99]
#     bounds = [(0, 1), (0, 1)]
#     result = minimize(objective_function, initial_guesses, args=(chi, N1, N2),bounds=bounds)
#     phi1A, phi1B = result.x
#     return phi1A, phi1B

# # Define the parameters
# N1 = 1  # Degree of polymerization for polymer 1
# N2 = 1  # Degree of polymerization for polymer 2
# chi_values = np.linspace(2, 5, 100)  # Range of chi values (starting from 0.5 to avoid trivial solutions)

# # Lists to store equilibrium compositions
# phi1_phase1 = []
# phi1_phase2 = []

# # Calculate equilibrium compositions for different chi values
# for chi in chi_values:
#     try:
#         phi1A, phi1B = solve_equilibrium_compositions(chi, N1, N2)
#         phi1_phase1.append(phi1A)
#         phi1_phase2.append(phi1B)
#     except ValueError:
#         # Handle cases where no solution is found
#         phi1_phase1.append(np.nan)
#         phi1_phase2.append(np.nan)

# # Plot the results
# plt.figure(figsize=(10, 6))
# plt.plot(phi1_phase1, chi_values, label='Phase 1 composition $\phi_1$', color='blue')
# plt.plot(phi1_phase2, chi_values, label='Phase 2 composition $\phi_2$', color='red')
# plt.xlabel('Composition $\phi_1$')
# plt.ylabel('$\chi$')
# plt.title('Equilibrium Compositions of Binary Polymer Mixtures')
# plt.legend()
# plt.grid(True)
# plt.show()

# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.optimize import minimize



# def flory_huggins_deriv(Nchi, N1_N2, phi1):
#     return  np.log(phi1) - N1_N2*np.log(1 - phi1) + 1 - N1_N2 + Nchi * (1 - 2 * phi1)

# def objective_function(vars, Nchi, N1_N2):
#     phi1A, phi1B = vars
#     muA = flory_huggins_deriv(Nchi, N1_N2, phi1A)
#     muB = flory_huggins_deriv(Nchi, N1_N2, phi1B)
#     return (muA - muB)**2

# def solve_equilibrium_compositions(Nchi, N1_N2, initial_guesses):
#     bounds = [(0, 1), (0, 1)]
#     result = minimize(objective_function, initial_guesses, args=(Nchi, N1_N2), bounds=bounds)
#     phi1A, phi1B = result.x
#     return phi1A, phi1B

# # Define the parameters
# N1_N2 = 1 
# nvals = 100
# chi_values = np.linspace(2, 6, nvals)  # Range of chi values (starting from 0.5 to avoid trivial solutions)

# # Lists to store equilibrium compositions
# phi1_phase1 = []
# phi1_phase2 = []

# # Calculate equilibrium compositions for different chi values
# initial_guesses = [0.49, 0.51]
# for Nchi in chi_values:
#     try:
#         phi1A, phi1B = solve_equilibrium_compositions(Nchi, N1_N2, initial_guesses)
#         phi1_phase1.append(phi1A)
#         phi1_phase2.append(phi1B)
        
#         # Update initial guesses for the next iteration
#         initial_guesses = [max(1e-8, phi1A - 0.001), min(1 - 1e-8, phi1B + 0.001)]
#     except ValueError:
#         # Handle cases where no solution is found
#         phi1_phase1.append(np.nan)
#         phi1_phase2.append(np.nan)

# # Plot the results
# plt.figure(figsize=(10, 6))
# plt.plot(phi1_phase1, chi_values, label='Phase 1 composition $\phi_1$', color='blue')
# plt.plot(phi1_phase2, chi_values, label='Phase 2 composition $\phi_2$', color='red')
# plt.xlabel('Composition $\phi_1$')
# plt.ylabel('$\chi$')
# plt.legend()
# plt.show()

def binodal(phi,chi,N1,N2):
    phiA = phi[0]
    phiB = phi[1]

    F = [0.,0.]

    F[0] = ((1/N1)*(1+np.log(phiA)) + chi*(1-phiA)) -  ((1/N1)*(1+np.log(phiB)) - chi*(1-phiB))
    F[1] = (1/N2)*(1+np.log(1-phiA)) +chi*(phiA) - (1/N2)*(1+np.log(1-phiB)) - chi*(phiA)

    return F

import numpy as np

from scipy.optimize import root

from matplotlib import pyplot as plt

import copy

 

def spinodal(chi,NA,NB):

    A = NB

    B = (NA-NB-2*chi*NA*NB)

    C = 2*chi*NA*NB

    return np.roots([C,B,A])

 

def binodal(phi,chi,NA,NB):

    phi1 = phi[0]

    phi2 = phi[1]

    F = [0.,0.]

    F[0] = -(phi1-phi2)*(1/NB-1/NA) + 1/NA*np.log((1-phi1)/(1-phi2))+chi*(phi1**2-phi2**2)

    F[1] = -(phi1-phi2)*(1/NB-1/NA) + 1/NB*np.log(phi1/phi2)+chi*((1-phi1)**2-(1-phi2)**2)

    return F

def fh_deriv(phi,chi,N1,N2):
    df = (1/N1)*np.log(phi) + (1/N1) - (1/N2)*np.log(1-phi) - (1/N2) -2*chi*phi + chi
    return df

def osmotic(phi,chi,N1,N2):
    osmo = -phi*((1/N1)-(1/N2)) -(1/N2)*np.log(1-phi) - chi*phi**2
    return osmo

# def binodal(phi,chi,N1,N2):
#     phiA = phi[0]
#     phiB = phi[1]

#     F = [0.,0.]

#     F[0] = fh_deriv(phiA,chi,N1,N2) - fh_deriv(phiB,chi,N1,N2)
#     F[1] = osmotic(phiA,chi,N1,N2) - osmotic(phiB,chi,N1,N2)

#     return F

 

NA = 1

NB = 2

 

phic = np.sqrt(NB)/(np.sqrt(NA)+np.sqrt(NB))

chic = (np.sqrt(NA)+np.sqrt(NB))**2/(2*NA*NB)

 

chi = np.linspace(20.,1.,50000)

phi1bn = np.zeros(len(chi))

phi2bn = np.zeros(len(chi))

phi1sp = np.zeros(len(chi))

phi2sp = np.zeros(len(chi))

 

x0 = [1e-4,1-1e-4]

 

for i in range(len(chi)):

   

    phisp = spinodal(chi[i]*chic,NA,NB)

    phi1sp[i] = phisp[0]

    phi2sp[i] = phisp[1]

 

    vals = root(lambda x: binodal(x,chi[i]*chic,NA,NB),x0)

    phi1bn[i] = vals.x[0]

    phi2bn[i] = vals.x[1]

    x0 = copy.copy(vals.x)

 

plt.plot(1-phi1bn,chi**-1,'b',label='Binodal')

plt.plot(1-phi2bn,chi**-1,'b',label='')

plt.plot(phi1sp,chi**-1,'r',label='Spinodal')

plt.plot(phi2sp,chi**-1,'r',label='')

plt.plot(phic,[1.],'ko',label='Critical')

plt.xlabel(r'$\phi$',fontsize=16)

plt.ylabel(r'$\chi_c/\chi$',fontsize=16)

plt.xlim(0,1)

plt.ylim(0.0,1.1)

plt.legend(frameon=False)

plt.show()
