import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)


#Read csv
df_99 = pd.read_csv("./Solidfrac_99.csv")
df_100 = pd.read_csv("./Solidfrac_100.csv")
df_101 = pd.read_csv("./Solidfrac_101.csv")

fig1,ax1 = plt.subplots(1,1,num=1)
ax1.plot(df_99["Time"],df_99["avg(SolidFrac)"],label=r"$r_{0} = 0.99 r^{*}$")
ax1.plot(df_100["Time"],df_100["avg(SolidFrac)"],label=r"$r_{0} = r^{*}$")
ax1.plot(df_101["Time"],df_101["avg(Solid_Frac)"],label=r"$r_{0} = 1.01 r^{*}$")
ax1.set_ylabel(r"$r$")
ax1.set_xlabel("Time")
ax1.legend()
plt.tight_layout()
plt.savefig("r.png",dpi=300)

plt.show()