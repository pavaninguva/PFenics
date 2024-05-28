import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import fsolve


#formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)

#Domain length
L = 5.0
nx = 101
dx = L/(nx-1)
xvals = np.linspace(0,L,nx)

#Timestepping parameters
tf = 10
nsteps = 10000
dt = tf/nsteps

#Material parameters
chi = 6.0
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


def backward_euler_fd_res(phi_new, phi_old):
    # Define discretized equations as F(phi_new) = 0
    res = np.zeros_like(phi_old)
    
    # Handle internal points
    for i in range(2, nx-2):
        mu_i__1_new = dfdphi(phi_new[i - 1]) - (kappa / (dx**2)) * (phi_new[i] - 2 * phi_new[i - 1] + phi_new[i - 2])
        mu_i_new = dfdphi(phi_new[i]) - (kappa / (dx**2)) * (phi_new[i + 1] - 2 * phi_new[i] + phi_new[i - 1])
        mu_i1_new = dfdphi(phi_new[i + 1]) - (kappa / (dx**2)) * (phi_new[i + 2] - 2 * phi_new[i + 1] + phi_new[i])
        res[i] = phi_new[i] - phi_old[i] - (dt / (dx**2)) * (M_func_half(phi_new[i], phi_new[i + 1]) * (mu_i1_new - mu_i_new) - M_func_half(phi_new[i], phi_new[i - 1]) * (mu_i_new - mu_i__1_new))

    # Handle Neumann boundary conditions
    # Left boundary (index 0 and 1)
    mu_0_new = dfdphi(phi_new[0]) - (2 * kappa / dx**2) * (phi_new[1] - phi_new[0])
    mu_1_new = dfdphi(phi_new[1]) - (kappa / dx**2) * (phi_new[2] - 2 * phi_new[1] + phi_new[0])
    res[0] = phi_new[0] - phi_old[0] - (2 * dt / dx**2) * (M_func_half(phi_new[0], phi_new[1]) * (mu_1_new - mu_0_new))
    
    mu_2_new = dfdphi(phi_new[2]) - (kappa / dx**2) * (phi_new[3] - 2 * phi_new[2] + phi_new[1])
    res[1] = phi_new[1] - phi_old[1] - (dt / dx**2) * (M_func_half(phi_new[1], phi_new[2]) * (mu_2_new - mu_1_new) - M_func_half(phi_new[1], phi_new[0]) * (mu_1_new - mu_0_new))

    # Right boundary (index nx-2 and nx-1)
    mu_N_new = dfdphi(phi_new[-1]) - (2 * kappa / dx**2) * (phi_new[-2] - phi_new[-1])
    mu_N_1_new = dfdphi(phi_new[-2]) - (kappa / dx**2) * (phi_new[-3] - 2 * phi_new[-2] + phi_new[-1])
    res[-1] = phi_new[-1] - phi_old[-1] - (2 * dt / dx**2) * (M_func_half(phi_new[-1], phi_new[-2]) * (mu_N_new - mu_N_1_new))
    
    mu_N_2_new = dfdphi(phi_new[-3]) - (kappa / dx**2) * (phi_new[-4] - 2 * phi_new[-3] + phi_new[-2])
    res[-2] = phi_new[-2] - phi_old[-2] - (dt / dx**2) * (M_func_half(phi_new[-2], phi_new[-1]) * (mu_N_new - mu_N_1_new) - M_func_half(phi_new[-2], phi_new[-3]) * (mu_N_1_new - mu_N_2_new))

    return res



# Define initial condition as a sigmoid function
def sigmoid(x, a, b):
    return 1 / (1 + np.exp(-a * (x - b)))

# Parameters for the sigmoid function
a = 5  # Controls the steepness of the sigmoid
b = L / 2  # Center of the sigmoid function
phi0 = 0.2 + 0.6 * sigmoid(xvals, a, b)
# phi0 = np.random.normal(loc=0.5, scale=0.05, size=nx)
phi = phi0.copy()



# Animation setup
fig, ax = plt.subplots()
line, = ax.plot(xvals, phi, label=r"$\phi$")
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
ax.set_xlim(0.0, L)
ax.set_ylim(0.0, 1.0)
ax.set_xlabel("x")
ax.set_ylabel(r"$\phi$")
ax.legend()

# Main time-stepping loop
def update(n):
    global phi
    phi_old = phi.copy()

    def wrapped_residual(phi_new):
        return backward_euler_fd_res(phi_new, phi_old)
    
    phi_new = fsolve(wrapped_residual, phi_old)
    phi = phi_new.copy()
    
    line.set_ydata(phi)
    time_text.set_text(f'Time = {n * dt:.2f}')
    return line, time_text

ani = animation.FuncAnimation(fig, update, frames=nsteps, blit=True, repeat=False)

plt.show()







