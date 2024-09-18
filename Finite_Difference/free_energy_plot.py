import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)

def flory_huggins(phi,chi):
    f = phi*np.log(phi) + (1-phi)*np.log(1-phi) + phi*(1-phi)*chi
    return f

def dfdphi(phi,chi):
    return -2*chi*phi + chi - np.log(1-phi) + np.log(phi)

#Compute spline_pot
def spline_pot(chi):
    spline_pot_fun = interpolate.CubicSpline(np.linspace(0,1,200),dfdphi(np.linspace(1e-16,1-1e-16,200),chi))
    return spline_pot_fun



# def flory_huggins2(phi,chi):
#     f = (1+phi)*np.log(1+phi) + (1-phi)*np.log(1-phi) - chi*phi**2
#     return f

# def quartic_pot(phi,A):
#     f = A*(phi**2)*(1-phi)**2
#     return f

# ## Approx for log part
# def log_d(phi,delta = 1e-25):
#     def func_cond(x):
#         return x <= delta
#     def func1(x):
#         return np.log(delta) + (phi-delta)/delta - ((phi-delta)**2)/(2*delta**2) + ((phi-delta)**3)/(3*delta**3)
#     def func2(x):
#         return np.log(x)

#     f = np.where(func_cond(phi),func1(phi),func2(phi))
#     return f

# def dfdphi_approx(phi,chi):
#     f = phi*log_d(phi) + (1.0-phi)*log_d(1.0-phi) + chi*phi*(1.0-phi)
#     return f

phi_vals = np.linspace(1e-16,1.0-1e-16,201)
# phi_vals2 = np.linspace(-1.0+1e-12,1.0-1e-12,501)

# #Plot
# fig1,ax1 = plt.subplots(1,3,figsize=(15,5))
# #Classic FH
# ax1[0].plot(phi_vals,flory_huggins(phi_vals,1.5),label=r"$\chi = 1.5$")
# ax1[0].plot(phi_vals,flory_huggins(phi_vals,2.0),label=r"$\chi = 2.0$")
# ax1[0].plot(phi_vals,flory_huggins(phi_vals,3.0),label=r"$\chi = 3.0$")
# ax1[0].plot(phi_vals,flory_huggins(phi_vals,6.0),label=r"$\chi = 6.0$")
# ax1[0].set_ylabel(r"$f$")
# ax1[0].set_xlabel(r"$\phi$")
# ax1[0].legend()
# ax1[0].set_title("Classic Flory-Huggins")
# #Scaled FH
# ax1[1].plot(phi_vals2,flory_huggins2(phi_vals2,1.0),label=r"$\chi = 0.5$")
# ax1[1].plot(phi_vals2,flory_huggins2(phi_vals2,1.5),label=r"$\chi = 1.5$")
# ax1[1].plot(phi_vals2,flory_huggins2(phi_vals2,2.0),label=r"$\chi = 2.0$")
# ax1[1].plot(phi_vals2,flory_huggins2(phi_vals2,3.0),label=r"$\chi = 3.0$")
# ax1[1].set_ylabel(r"$f$")
# ax1[1].set_xlabel(r"$\phi$")
# ax1[1].legend()
# ax1[1].set_title("Rescaled Flory-Huggins")
# #Quartic
# ax1[2].plot(phi_vals,quartic_pot(phi_vals,0.5),label=r"$A = 0.5$")
# ax1[2].plot(phi_vals,quartic_pot(phi_vals,1.0),label=r"$A = 1.0$")
# ax1[2].plot(phi_vals,quartic_pot(phi_vals,2.0),label=r"$A = 2.0$")
# ax1[2].plot(phi_vals,quartic_pot(phi_vals,5.0),label=r"$A = 5.0$")
# ax1[2].set_ylabel(r"$f$")
# ax1[2].set_xlabel(r"$\phi$")
# ax1[2].legend()
# ax1[2].set_title("Quartic Potential")

# fig1.tight_layout()

fig2,ax2 = plt.subplots(1,1,figsize=(6,6))
ax2.plot(phi_vals,dfdphi(phi_vals,3.0),"ok", label=r"FH, $\chi$ = 3.0",markerfacecolor='none')
spline_fun1 = spline_pot(3.0)
ax2.plot(phi_vals,spline_fun1(phi_vals),"-k",label=r"CS, $\chi$ = 3.0")

ax2.plot(phi_vals,dfdphi(phi_vals,50.0),"ob", label=r"FH, $\chi$ = 50.0",markerfacecolor='none')
spline_fun2 = spline_pot(50.0)
ax2.plot(phi_vals,spline_fun2(phi_vals),"-b",label=r"CS, $\chi$ = 50.0")

ax2.set_xlabel(r'$\phi$')
ax2.set_ylabel(r'$\frac{df}{d\phi}$')
ax2.legend()
fig2.tight_layout()

plt.savefig("Spline.png",dpi=300)

plt.show()
