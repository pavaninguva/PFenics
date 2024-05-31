import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import root
from scipy import interpolate


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
nsteps = 1000
dt = tf/nsteps

#Material parameters
chi = 50.0
kappa = (2/3)*chi

#Provide analytical expressions for f and dfdphi
def potential(phi):
    return phi*np.log(phi) + (1.0-phi)*np.log(1.0-phi) + chi*phi*(1.0-phi)

def dfdphi(phi):
    return -2*chi*phi + chi - np.log(1-phi) + np.log(phi)

def dfdphi_convex(phi,delta=1e-8):
    return - np.log(1-phi) + np.log(phi)

def dfdphi_concave(phi):
    return -2*chi*phi + chi

#Compute spline_pot
spline_pot = interpolate.CubicSpline(np.linspace(0,1,200),dfdphi(np.linspace(1e-16,1-1e-16,200)))



## Approx for log part
def log_d(phi,delta = 1e-8):
    # def func_cond(x):
    #     return np.abs(x) < delta  np.abs(x) > 1-delta
    def func1(x):
        return np.log(delta) + (phi-delta)/delta - ((phi-delta)**2)/(2*delta**2)
    def func2(x):
        return np.log(x)

    f = np.where((np.abs(phi)<delta)&(np.abs(phi) > 1-delta),func1(phi),func2(phi))
    return f

def dfdphi_approx(phi):
    f = -2*chi*phi + chi - log_d(1-phi) + log_d(phi)
    return f

def dfdphi_convex_approx(phi):
    return -log_d(1-phi) + log_d(phi)


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


def backward_euler_fd_res(phi_new, phi_old):
    #Define chem_pot
    mu_new = np.zeros_like(phi_old)

    mu_new[0] = dfdphi(phi_new[0]) -(2*kappa/(dx**2))*(phi_new[1] - phi_new[0])
    mu_new[-1] = dfdphi(phi_new[-1]) - (2*kappa/(dx**2))*(phi_new[-2]-phi_new[-1])
    mu_new[1:-1] = dfdphi(phi_new[1:-1]) - (kappa/(dx**2))*(phi_new[2:] -2*phi_new[1:-1] + phi_new[:-2])
    
    # Define discretized equations as F(phi_new) = 0
    res = np.zeros_like(phi_old)

    res[0] = (phi_new[0] - phi_old[0])/dt - (2/(dx**2))*(M_func_half(phi_new[0],phi_new[1]))*(mu_new[1] - mu_new[0])
    res[-1] = (phi_new[-1] - phi_old[-1])/dt - (2/(dx**2))*(M_func_half(phi_new[-1],phi_new[-2]))*(mu_new[-2]-mu_new[-1])
    res[1:-1] = (phi_new[1:-1] - phi_old[1:-1])/dt - (1/(dx**2))*(M_func_half(phi_new[1:-1], phi_new[2:])*(mu_new[2:]-mu_new[1:-1]) - M_func_half(phi_new[1:-1],phi_new[:-2])*(mu_new[1:-1]-mu_new[:-2]))

    return res

def theta_fd_res(phi_new, phi_old, theta=0.5):
    #The theta scheme introduces a parameter theta which can be varied from 0 - 1 
    #Theta = 1: backward euler
    #Theta = 0: Forward euler
    #Theta = 0.5: Crank-Nicolson

    #Define chem_pot
    mu_new = np.zeros_like(phi_old)
    mu_old = np.zeros_like(phi_old)

    mu_new[0] = dfdphi(phi_new[0]) -(2*kappa/(dx**2))*(phi_new[1] - phi_new[0])
    mu_new[-1] = dfdphi(phi_new[-1]) - (2*kappa/(dx**2))*(phi_new[-2]-phi_new[-1])
    mu_new[1:-1] = dfdphi(phi_new[1:-1]) - (kappa/(dx**2))*(phi_new[2:] -2*phi_new[1:-1] + phi_new[:-2])

    mu_old[0] = dfdphi(phi_old[0]) -(2*kappa/(dx**2))*(phi_old[1] - phi_old[0])
    mu_old[-1] = dfdphi(phi_old[-1]) - (2*kappa/(dx**2))*(phi_old[-2]-phi_old[-1])
    mu_old[1:-1] = dfdphi(phi_old[1:-1]) - (kappa/(dx**2))*(phi_old[2:] -2*phi_old[1:-1] + phi_old[:-2])

    res = np.zeros_like(phi_old)

    res[0] =  (phi_new[0] - phi_old[0])/dt - (theta*((2/(dx**2))*(M_func_half(phi_new[0],phi_new[1]))*(mu_new[1] - mu_new[0])) +
                                              (1-theta)*(2/(dx**2))*(M_func_half(phi_old[0],phi_old[1]))*(mu_old[1] - mu_old[0])
                                              )
    res[-1] = (phi_new[-1] - phi_old[-1])/dt -(theta*((2/(dx**2))*(M_func_half(phi_new[-1],phi_new[-2]))*(mu_new[-2]-mu_new[-1])) +
                                               (1-theta)*((2/(dx**2))*(M_func_half(phi_old[-1],phi_old[-2]))*(mu_old[-2]-mu_old[-1]))
                                                )
    res[1:-1] = (phi_new[1:-1] - phi_old[1:-1])/dt - (theta*((1/(dx**2))*(M_func_half(phi_new[1:-1], phi_new[2:])*(mu_new[2:]-mu_new[1:-1]) - M_func_half(phi_new[1:-1],phi_new[:-2])*(mu_new[1:-1]-mu_new[:-2]))) +
                                                      (1-theta)*((1/(dx**2))*(M_func_half(phi_old[1:-1], phi_old[2:])*(mu_old[2:]-mu_old[1:-1]) - M_func_half(phi_old[1:-1],phi_old[:-2])*(mu_old[1:-1]-mu_old[:-2])))
                                                        )
    
    return res


def convex_split_ch(phi_new, phi_old):
    #Define chem_pot
    mu_new = np.zeros_like(phi_old)

    mu_new[0] = dfdphi_convex(phi_new[0]) + dfdphi_concave(phi_old[0]) -(2*kappa/(dx**2))*(phi_new[1] - phi_new[0])
    mu_new[-1] = dfdphi_convex(phi_new[-1]) + dfdphi_concave(phi_old[-1]) - (2*kappa/(dx**2))*(phi_new[-2]-phi_new[-1])
    mu_new[1:-1] = dfdphi_convex(phi_new[1:-1]) + dfdphi_concave(phi_old[1:-1]) - (kappa/(dx**2))*(phi_new[2:] -2*phi_new[1:-1] + phi_new[:-2])
    
    # Define discretized equations as F(phi_new) = 0
    res = np.zeros_like(phi_old)

    res[0] = (phi_new[0] - phi_old[0])/dt - (2/(dx**2))*(M_func_half(phi_new[0],phi_new[1]))*(mu_new[1] - mu_new[0])
    res[-1] = (phi_new[-1] - phi_old[-1])/dt - (2/(dx**2))*(M_func_half(phi_new[-1],phi_new[-2]))*(mu_new[-2]-mu_new[-1])
    res[1:-1] = (phi_new[1:-1] - phi_old[1:-1])/dt - (1/(dx**2))*(M_func_half(phi_new[1:-1], phi_new[2:])*(mu_new[2:]-mu_new[1:-1]) - M_func_half(phi_new[1:-1],phi_new[:-2])*(mu_new[1:-1]-mu_new[:-2]))

    return res

def backward_euler_log_approx_res(phi_new, phi_old):
    #Define chem_pot
    mu_new = np.zeros_like(phi_old)

    mu_new[0] = dfdphi_approx(phi_new[0]) -(2*kappa/(dx**2))*(phi_new[1] - phi_new[0])
    mu_new[-1] = dfdphi_approx(phi_new[-1]) - (2*kappa/(dx**2))*(phi_new[-2]-phi_new[-1])
    mu_new[1:-1] = dfdphi_approx(phi_new[1:-1]) - (kappa/(dx**2))*(phi_new[2:] -2*phi_new[1:-1] + phi_new[:-2])
    
    # Define discretized equations as F(phi_new) = 0
    res = np.zeros_like(phi_old)

    res[0] = (phi_new[0] - phi_old[0])/dt - (2/(dx**2))*(M_func_half(phi_new[0],phi_new[1]))*(mu_new[1] - mu_new[0])
    res[-1] = (phi_new[-1] - phi_old[-1])/dt - (2/(dx**2))*(M_func_half(phi_new[-1],phi_new[-2]))*(mu_new[-2]-mu_new[-1])
    res[1:-1] = (phi_new[1:-1] - phi_old[1:-1])/dt - (1/(dx**2))*(M_func_half(phi_new[1:-1], phi_new[2:])*(mu_new[2:]-mu_new[1:-1]) - M_func_half(phi_new[1:-1],phi_new[:-2])*(mu_new[1:-1]-mu_new[:-2]))

    return res

def convexity_split_log_approx_res(phi_new, phi_old):
#     #Define chem_pot
    mu_new = np.zeros_like(phi_old)

    mu_new[0] = dfdphi_convex_approx(phi_new[0]) + dfdphi_concave(phi_old[0]) -(2*kappa/(dx**2))*(phi_new[1] - phi_new[0])
    mu_new[-1] = dfdphi_convex_approx(phi_new[-1]) + dfdphi_concave(phi_old[-1]) - (2*kappa/(dx**2))*(phi_new[-2]-phi_new[-1])
    mu_new[1:-1] = dfdphi_convex_approx(phi_new[1:-1]) + dfdphi_concave(phi_old[1:-1]) - (kappa/(dx**2))*(phi_new[2:] -2*phi_new[1:-1] + phi_new[:-2])
    
    # Define discretized equations as F(phi_new) = 0
    res = np.zeros_like(phi_old)

    res[0] = (phi_new[0] - phi_old[0])/dt - (2/(dx**2))*(M_func_half(phi_old[0],phi_old[1]))*(mu_new[1] - mu_new[0])
    res[-1] = (phi_new[-1] - phi_old[-1])/dt - (2/(dx**2))*(M_func_half(phi_old[-1],phi_old[-2]))*(mu_new[-2]-mu_new[-1])
    res[1:-1] = (phi_new[1:-1] - phi_old[1:-1])/dt - (1/(dx**2))*(M_func_half(phi_old[1:-1], phi_old[2:])*(mu_new[2:]-mu_new[1:-1]) - M_func_half(phi_old[1:-1],phi_old[:-2])*(mu_new[1:-1]-mu_new[:-2]))

    return res

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

# Parameters for the sigmoid function
a = 5  # Controls the steepness of the sigmoid
b = L / 2  # Center of the sigmoid function
phi0 = 0.2 + 0.6* sigmoid(xvals, a, b)
# phi0 = np.random.normal(loc=0.5, scale=0.05, size=nx)
phi = phi0.copy()


# Animation setup
fig, ax = plt.subplots()
line, = ax.plot(xvals, phi, label=r"$\phi$")
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
ax.set_xlim(0.0, L)
ax.set_ylim(-0.2, 1.2)
ax.set_xlabel("x")
ax.set_ylabel(r"$\phi$")
ax.hlines(y=0.0, xmin=0, xmax=5,color="r")
ax.hlines(y=1.0, xmin=0, xmax=5,color="r")
ax.legend()

# Main time-stepping loop
def update(n):
    global phi
    phi_old = phi.copy()
    # phi_old = np.clip(phi_old_,0+1e-8,1-1e-8)

    def wrapped_residual(phi_new):
        return backward_euler_splinepot_res(phi_new, phi_old)
    
    sol = root(wrapped_residual, phi_old,method="lm",tol=1e-10)
    phi_new = sol.x
    phi = phi_new.copy()
    # print(np.log(phi))
    
    line.set_ydata(phi)
    time_text.set_text(f'Time = {n * dt:.2f}')
    return line, time_text

ani = animation.FuncAnimation(fig, update, frames=nsteps, blit=True, repeat=False)

plt.show()







