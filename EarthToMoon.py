#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 10:27:11 2021

@author: whiffen
"""

import numpy as np

import matplotlib.pyplot as plt

class physicalQuantity:

    def __init__(self, value, error):
        self.value = value
        self.error = error



G = 6.7 * 10 ** (-11)  # N m2 /kg2

radius_earth = physicalQuantity(6356.752 * 10 ** 3, 0.0005 * 10 ** 3)  # m polar radius
mass_earth = physicalQuantity(5.9724 * 10 ** 24, 0.00005 * 10 ** 24)  # kg

radius_moon = physicalQuantity(1736.0 * 10 ** 3, 
        #0.05 * 10 ** 3
        5*10**3)  # m polar radius
mass_moon = physicalQuantity(0.07346 * 10 ** 24, 0.000005 * 10 * 24)  # kg

perigee = physicalQuantity(0.3633 * 10 ** 9, 0.00005 * 10 ** 9)  # m
apogee = physicalQuantity(0.4055 * 10 ** 9, 0.00005 * 10 ** 9)  # m

distance_EarthToMoon = 4 * 10 ** 8  # m

dt = 10# sec

v_0 = [3 * 10 **4, 0, 0]  # m/s
x_0 = [radius_earth.value + 100 * 10 ** 3, 0, 0]  # m 100km represents the Karman line

v_0 = np.array(v_0)
x_0 = np.array(x_0)  # m


class planet:

    def __init__(self, pos, radius, mass):
        self.pos = pos  # m
        self.radius = radius  # m
        self.mass = mass  # kg


class rocket:

    def __init__(self, pos, v, bodies):
        self.pos = pos
        self.v = v
        self.bodies = bodies

    def acceleration(self):
        a = 0
        for body in self.bodies:
            r = body.pos.value - self.pos
            r_abs = np.linalg.norm(r)
            a = a + ((G * body.mass.value) / (r_abs) ** 3) * r
        return a

    def move(self):
        self.a = self.acceleration()
        self.v = self.v + self.a * dt
        self.pos = self.pos + self.v * dt


# Let's calculate 

pos_earth = np.array([0, 0, 0], dtype='f')
# dtype = float32 so 24 bits are
# for mantissa and we are correct to at least 7 decimal places.
pos_earth = physicalQuantity(pos_earth, np.array([1, 1, 1]) * 0.5 * 10 ** -7)  # m
earth = planet(pos_earth, radius_earth, mass_earth)

pos_moon = [0.5 * (perigee.value + apogee.value), 0, 0]
pos_moon = physicalQuantity(np.array(pos_moon), np.array([1, 1, 1]) * (perigee.error + apogee.error))
moon = planet(pos_moon, radius_moon, mass_moon)


def hit_detected(object_1_pos, object_2_pos, d, time_taken):
    if time_taken == 0:
        if np.linalg.norm(object_1_pos - object_2_pos) < d:
            return True
        else:
            return False
    else:
        return False


def travel_to_the_moon(x_0, v_0):
    bodies = [earth, moon]
    my_rocket = rocket(x_0, v_0, bodies)
    t = 0
    time_taken = np.zeros((2,2))
    foo = 0

    while (t < 1 * 10 ** 8):
        distance_rocketToEarth= np.linalg.norm(my_rocket.pos - earth.pos.value)
        if distance_rocketToEarth > 2 * apogee.value:
            #print("Out of the vicinity of the Earth-Moon system")
            return -1

        if (distance_rocketToEarth < earth.radius.value):
            #print("Fell back to Earth!")
            return -2

        moon_pos = np.array([moon.pos.value - moon.pos.error, moon.pos.value + moon.pos.error ])
        d = np.array([radius_moon.value + radius_moon.error, radius_moon.value - radius_moon.error])


        for i  in range(0, 2): 
            for j in range(0,2):
                if (hit_detected(my_rocket.pos, moon_pos[i], d[j], time_taken[i][j])):
                    time_taken[i][j] = t
                    foo = foo + 1
                    
        if foo > 3:
            return time_taken

        # if (np.linalg.norm(my_rocket.pos - (moon.pos.value-moon.pos.error)) < moon.radius.value+moon.radius.error):
        # print("We have arrived to the moon! \nIt took us", round(t,2), "seconds or ", round(t/3600,3) , "hours.")
        # print("Initial speed was ", round(np.linalg.norm(v_0),2), "m/s or ", round(np.linalg.norm(v_0)*(3600/1000), 2), "km/h")
        #   time_taken_1 = t

        # if (np.linalg.norm(my_rocket.pos - (moon.pos.value-moon.pos.error)) < moon.radius.value-moon.radius.error):
        #   time_taken_2 = t

        else:
            #if (t%100==0):
                #print("The time",t,", the position", my_rocket.pos)
            my_rocket.move()
            t = t + dt
    return -99



def timeEstimation(x_0,v_0):

    
    time_taken = travel_to_the_moon(x_0,v_0)

    if isinstance(time_taken,int):
        if time_taken == -1:
             return -1
        if time_taken == -2:
             return -2
        else:
             return -99

    short_time = np.min(time_taken)
    long_time = np.max(time_taken)

    average_time = physicalQuantity( ((long_time+ 0.5*dt) +  (short_time - 0.5*dt))/2, long_time+0.5*dt-(short_time-0.5*dt) )

    return average_time
    



times = np.array([])
v_0_speeds = np.arange(12*10**3, 99*10**3, 10**3 )
e_1 = np.array([1,0,0]).reshape(3,1)
v_0_array = e_1 * v_0_speeds

for column in v_0_array.T:
    time_estimation = timeEstimation(x_0,column) 
    
    if isinstance(time_estimation, int):
        if time_estimation == -1:
            print("Rocket has left the vicinity of the Earth-Moon system with an initial velocity of", column)
        if time_estimation == -2:
            print("One of the initial velocities ", column ,"has a magnitude below the minimum initial speed", 11666)
            
        if time_estimation == -99:
            print("Unknown error has occurred with an initial velocity of ", column) 
    else:
        times = np.append(times,time_estimation)

avg_times = np.zeros(len(times))
error_times = np.zeros(len(times))

print(avg_times)

for i in range(len(times)):
    avg_times[i] = times[i].value
    error_times[i] = times[i].error

plt.plot(v_0_speeds, avg_times, 'r-')
plt.axis([12*10**3, 99*10**3, np.min(avg_times), np.max(avg_times)])
plt.show()


#b = timeEstimation(x_0,v_0)

#print("We have arrived to the moon! \nIt took us", round(b.value,2), "seconds or ", round(b.value/3600,3) , "hours.")





'''def main():
    
    t = 0

    while(True):
        
        if (np.linalg.norm(my_rocket.pos - earth.pos) < radius_earth):
            print("Didn't get to the moon!")
            break
        if (np.linalg.norm(my_rocket.pos - moon.pos) <= radius_moon):
            print("We have arrived to the moon! \nIt took us", round(t,2), "seconds or ", round(t/3600,2) , "hours.")
            print("Initial speed was ", round(v_0,2), "m/s or ", round(v_0*(3600/1000), 2), "km/h")
            break
        
        else:
            my_rocket.move()
            t = t + dt
    
    
    
if __name__ == "__main__": main()'''
