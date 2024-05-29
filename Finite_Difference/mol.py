import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp

#formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)

#Domain length
L = 5.0
nx = 51
dx = L/(nx-1)
xvals = np.linspace(0,L,nx)

#Timestepping parameters
tf = 2.0
nsteps = 10000
dt = tf/nsteps

#Material parameters
chi = 10.0
kappa = (2/3)*chi

#Provide analytical expressions for f and dfdphi
def potential(phi):
    return phi*np.log(phi) + (1.0-phi)*np.log(1.0-phi) + chi*phi*(1.0-phi)

def dfdphi(phi):
    return -2*chi*phi + chi - np.log(1-phi) + np.log(phi)

def mobility(phi):
    return phi*(1-phi)

def M_func_half(phi, phi_,option="1"):
    if option == "1":
        M_func = 0.5*(mobility(phi)+mobility(phi_))
    elif option == "2":
        M_func = mobility(0.5*(phi+phi_))
    elif option == "3":
        M_func = (2*mobility(phi)*mobility(phi_))/(mobility(phi) + mobility(phi_))
    else:
        raise ValueError("Specified option for mobility interpolation not available")
    return M_func


def CH_func(phi):
    #Define chem_pot
    mu = np.zeros_like(phi)

    mu[0] = dfdphi(phi[0]) -(2*kappa/(dx**2))*(phi[1] - phi[0])
    mu[-1] = dfdphi(phi[-1]) - (2*kappa/(dx**2))*(phi[-2]-phi[-1])
    mu[1:-1] = dfdphi(phi[1:-1]) - (kappa/(dx**2))*(phi[2:] -2*phi[1:-1] + phi[:-2])

    #define system of equations
    f = np.zeros_like(phi)

    f[0] = (2/(dx**2))*(M_func_half(phi[0],phi[1]))*(mu[1] - mu[0])
    f[-1] = (2/(dx**2))*(M_func_half(phi[-1],phi[-2]))*(mu[-2]-mu[-1])
    f[1:-1] = (1/(dx**2))*(M_func_half(phi[1:-1], phi[2:])*(mu[2:]-mu[1:-1]) - M_func_half(phi[1:-1],phi[:-2])*(mu[1:-1]-mu[:-2]))

    return f

#Wrap discretized system for ODE solver
def ode_system(t,u):
    print(t)
    return CH_func(u)

#Define initial condition:
# # Define initial condition as a sigmoid function
# def sigmoid(x, a, b):
#     return 1 / (1 + np.exp(-a * (x - b)))

# # Parameters for the sigmoid function
# a = 5  # Controls the steepness of the sigmoid
# b = L / 2  # Center of the sigmoid function
# phi0 = 0.2 + 0.6 * sigmoid(xvals, a, b)

phi0 = np.random.normal(loc=0.5, scale=0.05, size=nx)
#Solve 
sol = solve_ivp(ode_system, [0,tf], phi0, method="RK23", atol=1e-8, rtol=1e-6)

# Set up the figure and axis
fig, ax = plt.subplots()
line, = ax.plot(xvals, sol.y[:, 0], lw=2)
ax.set_xlim(0, L)
ax.set_ylim(-1, 1)
ax.set_xlabel('x')
ax.set_ylabel('u')
ax.set_xlim(0.0, L)
ax.set_ylim(-0.2, 1.2)
ax.hlines(y=0.0, xmin=0, xmax=5,color="r")
ax.hlines(y=1.0, xmin=0, xmax=5,color="r")
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

def animate(i):
    # line.set_ydata(sol.y[:, i])
    # time_text.set_text(f'Time = {sol.t[i]:.5f}')
    # return line, time_text
    idx = i * 50
    if idx < len(sol.t):
        line.set_ydata(sol.y[:, idx])
        time_text.set_text(f'Time = {sol.t[idx]:.5f}')
    return line, time_text

ani = animation.FuncAnimation(fig, animate, frames=len(sol.t), interval=0, blit=True)

plt.show()

