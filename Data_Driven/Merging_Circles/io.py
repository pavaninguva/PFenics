import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

plot_sample = True
sort_individual_files = False

# if sort_individual_files == True:
#     #Sort through each csv
#     for files in glob.glob("*.csv"):
#         pd_array = pd.read_csv("./%s" % files )
#         pd_array.sort_values(["Points:0","Points:1"],ascending=[False,True],inplace=True)
#         pd_array.drop(columns=["Points:2"],inplace=True)
#         counter_ = files.replace("t_","")
#         counter = round(int(counter_.replace(".csv",""))*0.2,1)
#         file_name = "t-"+str(counter)+".csv"
#         pd_array.to_csv(file_name,index=False)

if plot_sample == True:
    #Read in t=0
    t0 = pd.read_csv("./t-0.0.csv")
    #Plot t=0
    fig1 = plt.pcolor(t0["Points:0"].to_numpy().reshape(101,101),t0["Points:1"].to_numpy().reshape(101,101),t0["f_0"].to_numpy().reshape(101,101))

    plt.show()