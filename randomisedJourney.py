import numpy as np
import pandas as pd
import matplotlib as plt
from random import  uniform

G = 6.7 * 10 ** (-11) # universal gravitational constant 
dr = 3*10**5# (m) distance steps
h = 100*10**3 # (m) height of the atmosphere above the Earth's surface

class planet:

    def __init__(self, pos, radius, mass):
        self.pos = pos  # m
        self.radius = radius  # m
        self.mass = mass  # kg

#earth = planet(p_e, r_e, m_e)
#moon = planet(p_m, r_m, m_m)


class rocket:

    def __init__(self, pos, v, earth,moon):
    
        self.pos = pos
        self.v = v
        self.earth = earth
        self.moon = moon

    def acceleration(self,pos):
        a = 0
        bodies = [self.earth,self.moon]
        for body in bodies:
            r = body.pos - pos
            r_abs = abs(r)
            a = a + ((G * body.mass) / (r_abs) ** 3) * r
            
        return a

    def move(self):
        
        a_n = self.acceleration(self.pos)
        a_n1 = self.acceleration(self.pos + dr)

        self.v[0] = (self.v[0] + np.sqrt(self.v[0]**2 + 4*( (a_n ) )*dr))/2
        self.v[1] = (self.v[1] + np.sqrt(self.v[1]**2 + 4*( (a_n1) )*dr))/2
        self.pos = self.pos + dr 


    def getTimeIncrease(self):

        dt_small = (-0.1)*np.ones(2)
        dt_large = (-0.1)*np.ones(2)


        zeroAccelerationPoint = (self.moon.pos ) / ( np.sqrt(self.moon.mass / self.earth.mass) + 1 )

        for i in range(0,2):

            if self.pos < zeroAccelerationPoint and self.pos+dr < zeroAccelerationPoint:

                v_n = self.v[i]
                #a_n = self.a
                a_n1 = self.acceleration(self.pos+dr)
 
                dt_small[i] = dr / v_n
                dt_large[i] = ( -v_n + np.sqrt(v_n**2 + 4 * a_n1 * dr ) ) / ( 2 * a_n1 ) 

            if self.pos+dr == zeroAccelerationPoint:

                v_n = self.v[i]
                dt_small[i] =  dr / v_n
                dt_large[i] = dr / v_n

            if self.pos +dr > zeroAccelerationPoint and self.pos < zeroAccelerationPoint:

                x_n = self.pos
                v_n = self.v[i]
                a_n = self.acceleration(self.pos)
                a_n1 = self.acceleration(self.pos+dr)
            
                dt_small[i] = min( ( -v_n + np.sqrt(v_n**2 + 4 * a_n1 * dr ) ) / ( 2 * a_n1 ), dr / v_n) 
                dt_large[i] = dr / (v_n + a_n*(zeroAccelerationPoint -x_n)) 


            if self.pos + dr > zeroAccelerationPoint and self.pos >= zeroAccelerationPoint:

                v_n = self.v[i]
                a_n1 = self.acceleration(self.pos+dr)
                
                dt_small[i] =  ( -v_n + np.sqrt(v_n**2 + 4 * a_n1 * dr ) ) / ( 2 * a_n1 ) 
                dt_large[i] = dr / v_n

        return [min(dt_small), max(dt_large)]
            
            

def timeTaken(x_0,v_0,earth,moon):

    myrocket = rocket(x_0,v_0, earth,moon)

    t_smaller = 0
    t_larger = 0

    counter = 0

    while counter < 1*10**10:
        
        distancerocketearth = abs(myrocket.pos - earth.pos)
        distancerocketmoon = abs(myrocket.pos - moon.pos)

        if distancerocketearth > 2 * moon.pos:
            break
            return -99
        if distancerocketearth < earth.radius:
            break
            return -99
        if distancerocketmoon < moon.radius:
            break
            
        else:
            dt = myrocket.getTimeIncrease()
            if dt[0] < 0 or dt[1] < 0:
                print("going backwards in time")
                print((moon.pos-myrocket.pos)/moon.pos)
                break
            
            t_smaller = t_smaller + dt[0]
            t_larger = t_larger + dt[1]
            myrocket.move()
            counter = counter +1


    return [t_smaller, t_larger]

def randomise(variable, maxError):
    variable = uniform(variable+maxError, variable-maxError)
    return variable


def experimentTime(N, v_0):

    time = np.zeros((N,2))
    
    for j in range(N):
    
        p_m = randomise(0.4055 * 10 ** 9, 0.00005 * 10 ** 9)
        p_e = 0
    
        m_m = randomise(0.07346 * 10 ** 24, 0.000005 * 10 * 24)
        m_e = randomise(5.9724 * 10 ** 24, 0.00005 * 10 ** 24)
    
        r_m = randomise( 1736.0 * 10 ** 3, 0.05*10**3)
        r_e = randomise(6356.752 * 10 ** 3, 0.0005 * 10 ** 3)
    
        earth = planet(p_e, r_e, m_e)
        moon = planet(p_m, r_m, m_m)
    
        T = timeTaken(h + r_e, np.array([v_0, v_0]), earth, moon)
        time[j] = [T[0], T[1]]

    return [np.mean(time), np.std(time)]

#From the miminum_initial_speed.py file, we know that the minimum initial
# speed needed to reach the moon is between 10.9869 and 10.9871 km/s. Let's
# calculate the time it takes to reach the moon at 10 percent more than the 
# minimum initial speed.

rep = 20 # number of random events

v_min_low = 10.9869*10**3
v_min_high = 10.9871*10**3

v_low = 1.1*v_min_low
v_high = 1.1*v_min_high

T_10_percent = [experimentTime(rep,v_low), experimentTime(rep,v_high)]

print("The mean time to reach the moon was: ", T_10_percent[0][0]/3600, "and the standard error", T_10_percent[0][1]/3600, "hours for the lower estimate of 1.1*v_min")
print("The mean time to reach the moon was: ", T_10_percent[1][0]/3600, "and the standard error", T_10_percent[1][1]/3600, "hours for the higher  estimate of 1.1*v_min")
print("Let's calculate the time it takes to reach the moon")

# The following part is for calculating the time it takes to reach 
# the moon with many different initial speeds. This data will be used
# to create a graph of time against initial speeds, with error bars
# for the time. This graph can be found in the plot_time_vs_vel.py file.

bins = 10

v_0_array = np.linspace(12*10**3, 15*10**3, num=bins)

time_mean_array = np.zeros(bins)
time_sd_array   = np.zeros(bins)

rep =20

for i in range(bins):

    time = np.zeros((bins,2))

    time = experimentTime(rep, v_0_array[i])
    time_mean_array[i] = time[0]
    time_sd_array[i] = time[1]


# Let's save the data into a a separate file 

d = {'vel': v_0_array, 'mean time': time_mean_array, 'std time': time_sd_array}

time_vel_data= pd.DataFrame(d)
time_vel_data.to_csv("time_vel_data.csv")

