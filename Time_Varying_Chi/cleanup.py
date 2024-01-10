import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

plot_sample = False
sort_individual_files = True

if sort_individual_files == True:
    #Sort through each csv
    for files in glob.glob("*.csv"):
        pd_array = pd.read_csv("./%s" % files )
        pd_array.sort_values(["Points:0","Points:1"],ascending=[False,True],inplace=True)
        pd_array.drop(columns=["Points:2"],inplace=True)
        counter_ = files.replace("t_","")
        counter = round(int(counter_.replace(".csv",""))*2.0,1)
        file_name = "t-"+str(counter)+".csv"
        pd_array.to_csv(file_name,index=False)