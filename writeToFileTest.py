import numpy as np
import pandas as pd
import csv


time_mean_array = np.zeros(10)
time_sd_array = np.zeros(10)
time_mean_array[9] = 99

time = pd.DataFrame(time_mean_array, time_sd_array)

time.to_csv("time_values.csv")

'''
file = open('time_values.csv','w')

for i in range(10):
    file.write(time_mean_array[i])
    file.writ(time_sd_array[i])

file.close()
'''


