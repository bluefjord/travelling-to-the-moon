import numpy as np
import voyageToTheMoon as voyage
import matplotlib.pyplot as plt

h = 100*10**3
dt = 10

p_m = 0.4055 * 10 ** 9
p_e = 0 
m_m = 0.07346 * 10 ** 24
m_e = 5.9724 * 10 ** 24
r_m = 1736.0 * 10 ** 3
r_e = 6356.752 * 10 ** 3

earth = voyage.planet(p_e, r_e, m_e)
moon = voyage.planet(p_m, r_m, m_m)

bodies = [earth,moon]

timeIncrement = np.array([100,50,30,10,1])

N = 10

timeArray = np.zeros((len(timeIncrement), N))

initialSpeeds = np.linspace(12*10**3, 14*10**3, num=N)


i = 0
for dt in timeIncrement:
    j = 0
    for v in initialSpeeds:
        myrocket = voyage.rocket(h+r_e,v,bodies)
    
        time = voyage.timeTaken(myrocket,moon,earth,dt)
        timeArray[i][j] = time
        j = j + 1
    i = i +1

#print(timeArray)
print(timeArray[1])


fig = plt.figure()
for i in range(0,len(timeIncrement)):
    plt.plot(initialSpeeds/1000, timeArray[i]/3600, label='time increment = %d seconds' %timeIncrement[i] )
plt.title('Time taken vs initial speed')
plt.xlabel('Speed (km/s)')
plt.ylabel('Time (h)')
plt.legend()
plt.show()


