import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

time_vel_data = pd.read_csv('time_vel_data.csv')
print(time_vel_data.head()) # to display the first 5 lines of loaded data
 
time_vel_data = time_vel_data.to_numpy()
time_vel_data = np.transpose(time_vel_data)

v = time_vel_data[1]
t = time_vel_data[2]
t_error = time_vel_data[3]



#t_error = 1000*np.ones(len(v))


fig = plt.figure()
plt.errorbar(v*10**-3,t/3600,yerr=t_error/3600)
plt.title('Time taken vs speed')
plt.xlabel('Speed (km/s)')
plt.ylabel('Time (h)')
plt.show()
fig.savefig('graph_time_vs_vel.png')
