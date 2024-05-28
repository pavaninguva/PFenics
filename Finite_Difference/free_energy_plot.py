import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)

def flory_huggins(phi,chi):
    f = phi*np.log(phi) + (1-phi)*np.log(1-phi) + phi*(1-phi)*chi
    return f

def flory_huggins2(phi,chi):
    f = (1+phi)*np.log(1+phi) + (1-phi)*np.log(1-phi) - chi*phi**2
    return f

def quartic_pot(phi,A):
    f = A*(phi**2)*(1-phi)**2
    return f

phi_vals = np.linspace(1e-12,1.0-1e-12,501)
phi_vals2 = np.linspace(-1.0+1e-12,1.0-1e-12,501)

#Plot
fig1,ax1 = plt.subplots(1,3,figsize=(15,5))
#Classic FH
ax1[0].plot(phi_vals,flory_huggins(phi_vals,1.5),label=r"$\chi = 1.5$")
ax1[0].plot(phi_vals,flory_huggins(phi_vals,2.0),label=r"$\chi = 2.0$")
ax1[0].plot(phi_vals,flory_huggins(phi_vals,3.0),label=r"$\chi = 3.0$")
ax1[0].plot(phi_vals,flory_huggins(phi_vals,6.0),label=r"$\chi = 6.0$")
ax1[0].set_ylabel(r"$f$")
ax1[0].set_xlabel(r"$\phi$")
ax1[0].legend()
ax1[0].set_title("Classic Flory-Huggins")
#Scaled FH
ax1[1].plot(phi_vals2,flory_huggins2(phi_vals2,1.0),label=r"$\chi = 0.5$")
ax1[1].plot(phi_vals2,flory_huggins2(phi_vals2,1.5),label=r"$\chi = 1.5$")
ax1[1].plot(phi_vals2,flory_huggins2(phi_vals2,2.0),label=r"$\chi = 2.0$")
ax1[1].plot(phi_vals2,flory_huggins2(phi_vals2,3.0),label=r"$\chi = 3.0$")
ax1[1].set_ylabel(r"$f$")
ax1[1].set_xlabel(r"$\phi$")
ax1[1].legend()
ax1[1].set_title("Rescaled Flory-Huggins")
#Quartic
ax1[2].plot(phi_vals,quartic_pot(phi_vals,0.5),label=r"$A = 0.5$")
ax1[2].plot(phi_vals,quartic_pot(phi_vals,1.0),label=r"$A = 1.0$")
ax1[2].plot(phi_vals,quartic_pot(phi_vals,2.0),label=r"$A = 2.0$")
ax1[2].plot(phi_vals,quartic_pot(phi_vals,5.0),label=r"$A = 5.0$")
ax1[2].set_ylabel(r"$f$")
ax1[2].set_xlabel(r"$\phi$")
ax1[2].legend()
ax1[2].set_title("Quartic Potential")

fig1.tight_layout()

plt.show()
