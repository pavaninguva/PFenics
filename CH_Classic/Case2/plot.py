import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)

#Read CSVs
df1 = pd.read_csv("./CH2_100X100_dt002.csv")
df2 = pd.read_csv("./CH2_100X100_dt001.csv")
df3 = pd.read_csv("./CH2_200X200_dt001.csv")

fig1, ax1 = plt.subplots(1,1,num=1)
ax1.plot(df1["Time"],df1["c"], "-k")
ax1.plot(df2["Time"],df2["c"], "--k")
ax1.plot(df3["Time"],df3["c"], "-.k")
ax1.set_ylabel(r"$\bar{c}$",color="k")
ax1.tick_params(axis="y",color="k")
ax1.set_xlabel("Time")
ax1.set_ylim([-1.0,1.0])
ax11 = ax1.twinx()
ax11.plot(df1["Time"],df1["energy"],"-r",label=r"100$\times$100, dt = 0.02")
ax11.plot(df2["Time"],df2["energy"],"--r",label=r"100$\times$100, dt = 0.01")
ax11.plot(df3["Time"],df3["energy"],"-.r",label=r"200$\times$200, dt = 0.01")
ax11.set_ylabel(r"$F$ (Energy)",color="r")
ax11.tick_params(axis="y",labelcolor="r")
ax11.legend()
plt.tight_layout()
plt.savefig("energydecay.png",dpi=300)

plt.show()