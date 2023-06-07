import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)

#read CSV
df = pd.read_csv("./data.csv")

fig1,ax1 = plt.subplots(1,1,num=1)
ax1.plot(df["Time"],df["c"],color="k")
ax1.set_ylabel(r"$\bar{c}$",color="k")
ax1.set_xlabel("Time")
ax1.tick_params(axis="y",color="k")
ax1.set_ylim([0.0,1.0])
ax11 = ax1.twinx()
ax11.plot(df["Time"],df["energy"],color="r")
ax11.set_ylabel(r"$F$ (Energy)",color="r")
ax11.tick_params(axis="y",labelcolor="r")
plt.tight_layout()
plt.savefig("conservation.png",dpi=300)

plt.show()
