import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Read data.csv
df = pd.read_csv("./data.csv")

#formatting
plt.rcParams["text.usetex"] = True
plt.rc('font', family='serif')
plt.rc("font",size=14)

fig1,ax1 = plt.subplots(1,1,num=1)
ax1.plot(df["Time"],df["c"],color="k")
ax1.set_ylabel(r"$c$",color="k")
ax1.tick_params(axis="y",color="k")
ax1.set_xlabel("Time")
ax11 = ax1.twinx()
ax11.plot(df["Time"],df["energy"],color="r")
ax11.set_ylabel(r"$F$ (Energy)",color="r")
ax11.tick_params(axis="y",labelcolor="r")
plt.tight_layout()
plt.savefig("conservation.png",dpi=300)

#Read area.csv
# a_df = pd.read_csv("./area.csv")
# a200_df = pd.read_csv("./area200.csv")
# fig2,ax2 = plt.subplots(1,1,num=2)
# ax2.plot(a_df["Time"],a_df["area"],label=r"Simulation, $100\times 100$")
# ax2.plot(a200_df["Time"],a200_df["area"],label=r"Simulation, $200\times 200$")
# ax2.plot(a_df["Time"],a_df["ana"],label="Analytical")
# ax2.legend()
# ax2.set_ylabel("Area")
# ax2.set_xlabel("Time")
# plt.tight_layout()
# plt.savefig("area.png",dpi=300)

plt.show()
