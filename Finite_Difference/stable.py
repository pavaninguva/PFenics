import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import root
from scipy import interpolate
import os


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
tf = 1.0
nsteps = 500
dt = tf/nsteps

#Material parameters
chi = 50.0
kappa = (2/3)*chi

def dfdphi(phi):
    return -2*chi*phi + chi - np.log(1-phi) + np.log(phi)
#Compute spline_pot
spline_pot = interpolate.CubicSpline(np.linspace(0,1,200),dfdphi(np.linspace(1e-16,1-1e-16,200)))


def mobility(phi,option=1):
    if option == 1:
        mobility = phi*(1-phi)
    elif option ==2:
        mobility = 1
    return mobility

def M_func_half(phi, phi_,option="2"):
    if option == "1":
        M_func = 0.5*(mobility(phi)+mobility(phi_))
    elif option == "2":
        M_func = mobility(0.5*(phi+phi_))
    elif option == "3":
        M_func = (2*mobility(phi)*mobility(phi_))/(mobility(phi) + mobility(phi_))
    else:
        raise ValueError("Specified option for mobility interpolation not available")
    return M_func

def backward_euler_splinepot_res(phi_new, phi_old):
    #Define chem_pot
    mu_new = np.zeros_like(phi_old)

    mu_new[0] = spline_pot(phi_new[0]) -(2*kappa/(dx**2))*(phi_new[1] - phi_new[0])
    mu_new[-1] = spline_pot(phi_new[-1]) - (2*kappa/(dx**2))*(phi_new[-2]-phi_new[-1])
    mu_new[1:-1] = spline_pot(phi_new[1:-1]) - (kappa/(dx**2))*(phi_new[2:] -2*phi_new[1:-1] + phi_new[:-2])
    
    # Define discretized equations as F(phi_new) = 0
    res = np.zeros_like(phi_old)

    res[0] = (phi_new[0] - phi_old[0])/dt - (2/(dx**2))*(M_func_half(phi_new[0],phi_new[1]))*(mu_new[1] - mu_new[0])
    res[-1] = (phi_new[-1] - phi_old[-1])/dt - (2/(dx**2))*(M_func_half(phi_new[-1],phi_new[-2]))*(mu_new[-2]-mu_new[-1])
    res[1:-1] = (phi_new[1:-1] - phi_old[1:-1])/dt - (1/(dx**2))*(M_func_half(phi_new[1:-1], phi_new[2:])*(mu_new[2:]-mu_new[1:-1]) - M_func_half(phi_new[1:-1],phi_new[:-2])*(mu_new[1:-1]-mu_new[:-2]))

    return res

# Define initial condition as a sigmoid function
def sigmoid(x, a, b):
    return 1 / (1 + np.exp(-a * (x - b)))

# # Parameters for the sigmoid function
a = 5  # Controls the steepness of the sigmoid
b = L / 2  # Center of the sigmoid function
phi0 = 0.2 + 0.6* sigmoid(xvals, a, b)
# phi0 = np.random.normal(loc=0.5, scale=0.05, size=nx)
phi = phi0.copy()


# Animation setup
fig, (ax1, ax2,ax3,ax4) = plt.subplots(4, 1, figsize=(6, 10))
line1, = ax1.plot(xvals, phi, label=r"$\phi$")
time_text1 = ax1.text(0.05, 0.9, '', transform=ax1.transAxes)
ax1.set_xlim(0.0, L)
ax1.set_ylim(-0.2, 1.2)
ax1.set_xlabel("x")
ax1.set_ylabel(r"$\phi$")
ax1.hlines(y=0.0, xmin=0, xmax=5, color="r")
ax1.hlines(y=1.0, xmin=0, xmax=5, color="r")
ax1.legend()

integrated_phi = np.trapz(phi, xvals)
integrated_values = [integrated_phi]
line2, = ax2.plot([0], [integrated_phi], label=r"$Total Mass$")
time_text2 = ax2.text(0.05, 0.9, '', transform=ax2.transAxes)
ax2.set_xlim(0, tf)
ax2.set_ylim(0, np.max(integrated_values) * 1.1)
ax2.set_xlabel("Time")
ax2.set_ylabel("Total Mass")

log_approx = interpolate.CubicSpline(np.linspace(0,1,200),np.log(np.linspace(1e-18,1-1e-18,200)))


def calculate_integral(phi):
    term1 = phi * log_approx(phi)
    term2 = (1 - phi) * log_approx(1 - phi)
    term3 = chi * (1 - phi) * phi
    dphi_dx = np.gradient(phi, dx)
    term4 = (kappa / 2) * dphi_dx**2
    return np.trapz(term1 + term2 + term3 + term4, xvals)

integral_I = calculate_integral(phi)
integral_values = [integral_I]
line3, = ax3.plot([0], [integral_I], label=r"$Total Energy$")
time_text3 = ax3.text(0.05, 0.9, '', transform=ax3.transAxes)
ax3.set_xlim(0, tf)
ax3.set_ylim(0, np.max(integral_values) * 1.1)
ax3.set_xlabel("Time")
ax3.set_ylabel("Total Energy")


def surf_integral(phi):
    dphi_dx = np.gradient(phi,dx)
    term1 = kappa*dphi_dx**2
    return np.trapz(term1,xvals)

integral_surf = surf_integral(phi)
integral_surf_values = [integral_surf]
line4, = ax4.plot([0], [integral_surf], label=r"$Total Energy$")
time_text4 = ax4.text(0.05, 0.9, '', transform=ax4.transAxes)
ax4.set_xlim(0, tf)
ax4.set_ylim(0, np.max(integral_values) * 1.1)
ax4.set_xlabel("Time")
ax4.set_ylabel(r"$\sigma$")


fig.tight_layout()


# Create a directory for snapshots
snapshot_dir = 'snapshots'
os.makedirs(snapshot_dir, exist_ok=True)

# Snapshot interval
snapshot_interval = 25

# Main time-stepping loop
def update(n):
    global phi
    phi_old = phi.copy()

    def wrapped_residual(phi_new):
        return backward_euler_splinepot_res(phi_new, phi_old)
    
    sol = root(wrapped_residual, phi_old, method="lm", tol=1e-10)
    phi_new = sol.x
    phi = phi_new.copy()

    # Update the phi plot
    line1.set_ydata(phi)
    time_text1.set_text(f'Time = {n * dt:.2f}')

    # Integrate phi over the domain
    integrated_phi = np.trapz(phi, xvals)
    integrated_values.append(integrated_phi)

    # Update the integrated phi plot
    line2.set_data(np.linspace(0, (n + 1) * dt, len(integrated_values)), integrated_values)
    ax2.set_ylim(2.0, np.max(integrated_values) * 1.1)
    time_text2.set_text(f'Integrated Value = {integrated_phi:.8f}')

    # Calculate the additional integral I
    integral_I = calculate_integral(phi)
    integral_values.append(integral_I)

    # Update the additional integral plot
    line3.set_data(np.linspace(0, (n + 1) * dt, len(integral_values)), integral_values)
    ax3.set_ylim(0, np.max(integral_values) * 1.1)
    time_text3.set_text(f'Integral I = {integral_I:.5f}')

    # Calculate the additional integral I
    integral_surf = surf_integral(phi)
    integral_surf_values.append(integral_surf)

    # Update the additional integral plot
    line4.set_data(np.linspace(0, (n + 1) * dt, len(integral_surf_values)), integral_surf_values)
    ax4.set_ylim(0, np.max(integral_surf_values) * 1.1)
    time_text4.set_text(f'$\sigma = {integral_surf:.5f}$')

    # Save snapshot every nth frame
    if n % snapshot_interval == 0:
        plt.savefig(os.path.join(snapshot_dir, f'snapshot_{int(n/snapshot_interval)}.png'))

    return line1, time_text1, line2, time_text2, line3, time_text3, line4, time_text4

ani = animation.FuncAnimation(fig, update, frames=nsteps, blit=True, repeat=False)

plt.show()