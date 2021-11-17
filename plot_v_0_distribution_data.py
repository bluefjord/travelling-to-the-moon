import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

v_0_distribution_data = pd.read_csv('v_0_distribution_data.csv')

print(v_0_distribution_data)

v_0_distribution_data = v_0_distribution_data.to_numpy()

v_0_distribution_data = np.transpose(v_0_distribution_data)

v = v_0_distribution_data[1]
count = v_0_distribution_data[2]

fig = plt.figure()
plt.plot(v,count)
plt.title('Histogram of minimum initial velocities')
plt.xlabel('Minimum initial speeds (m/s) ')
plt.ylabel('Number of occurences')
plt.show()
fig.savefig('graph_v_0_distribution.png')

