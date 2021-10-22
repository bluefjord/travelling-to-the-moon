#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 10:27:11 2021

@author: whiffen
"""

import numpy as np


class physical_quantity:
    
    def __init__(self,value,error):
        self.value = value
        self.error = error


G = 6.7 * 10**(-11) # N m2 /kg2

radius_earth =  physical_quantity(6356.752 *10**3, 0.0005 *10**3) # m polar radius
mass_earth = physical_quantity(5.9724 *10**24, 0.00005 *10**24) # kg

radius_moon = physical_quantity(1736.0 *10**3, 0.05 *10**3) # m polar radius
mass_moon = physical_quantity(0.07346 *10**24, 0.000005 *10*24) # kg
                              
perigree = physical_quantity(0.3633 *10**9,0.00005 *10**9) # m
apogee = physical_quantity(0.4055 *10**9, 0.00005 *10**9) # m
 
                              
distance_EarthToMoon = 4 *10**8 # m

dt = 100 # sec

v_0 = [7*10**8,0,0]  # m/s
x_0 = [100*10*3,0,0] # m 100km represents the Karman line

v_0 = np.array(v_0)
x_0 = physical_quantity(np.array(x_0),0.5*10**-7) # m



class planet:
    
    def __init__(self,pos,radius,mass):
        self.pos = pos # m
        self.radius = radius# m
        self.mass = mass # kg


class rocket:
    
    def __init__(self, pos,v,bodies):
        self.pos = pos
        self.v = v
        self.bodies = bodies

    def acceleration(self):
        a = 0
        for body in self.bodies:
            r = body.pos - self.pos
            r_abs = np.linalg.norm(r)
            a = a + ( (G * body.mass)/(r_abs)**3 )*r
        return a
    
    def move(self):
        self.a = self.acceleration()
        self.v = self.v + self.a*dt
        self.pos = self.pos + self.v*dt
    
# Let's calculate 
        
pos_earth = np.array([0,0,0],dtype = 'f') # dtype = float32 so 24 bits are 
# for mantissa and we are correct to at least 7 decimal places.
pos_earth = physical_quantity(pos_earth, np.array([1,1,1])*0.5*10**-7 ) # m
earth = planet(pos_earth, radius_earth, mass_earth )

pos_moon = [0.5*(perigree.value+apogee.value),0,0]
pos_moon= physical_quantity(np.array(pos_moon), np.array([1,1,1])*(perigree.error + apogee.error))
moon = planet(pos_moon, radius_moon, mass_moon)


    
    
def hit_detected(object_1_pos, object_2_pos, d):
    if( np.linalg.norm(object_1_pos - object_2_pos) < d):
        return True
    else:
        return False
    

def travel_to_the_moon(x_0,v_0):
    bodies = [earth,moon]
    my_rocket = rocket(x_0, v_0, bodies)  
    t=0
    time_taken = np.ones((2,2))
    
    while(t < 1*10**8):
        distance_rocketToEarthAtmos = np.linalg.norm(my_rocket.pos.value - np.array([100*10**3,0,0]))
        if(distance_rocketToEarthAtmos > 2*apogee.value):
            print("Out of the vicinity of the Earth-Moon system")
            return -1
        
        if (distance_rocketToEarthAtmos < 0):
            print("Fell back to Earth!")
            return -1
        
        vec1 = np.array([my_rocket.pos.value+my_rocket.pos.error, my_rocket.pos.value-my_rocket.pos.error])
        vec2 = np.array([moon.pos.value + moon.pos.error, moon.pos.value-moon.pos.error ])
        d = np.array([radius_moon.value+radius_moon.error,radius_moon.value - radius_moon.error])
        
        foo = 0
        for i in range(0,1):
            for j in range(0,1):
                if(hit_detected(vec1[i],vec2[j],d[0])):
                    foo = foo +1
                    time_taken[i][j] = t
                    
        if foo> 0:
            return time_taken
        
        #if (np.linalg.norm(my_rocket.pos - (moon.pos.value-moon.pos.error)) < moon.radius.value+moon.radius.error):
            #print("We have arrived to the moon! \nIt took us", round(t,2), "seconds or ", round(t/3600,3) , "hours.")
            #print("Initial speed was ", round(np.linalg.norm(v_0),2), "m/s or ", round(np.linalg.norm(v_0)*(3600/1000), 2), "km/h")
         #   time_taken_1 = t
            
        #if (np.linalg.norm(my_rocket.pos - (moon.pos.value-moon.pos.error)) < moon.radius.value-moon.radius.error):
         #   time_taken_2 = t
        
            
        else:
            my_rocket.move()
            t = t + dt
    return -99

travel_to_the_moon(x_0,v_0)
'''
# Let's find the minimum initial velocity needed to reach the moon
        
a = [1*10**4,0,0]# m/s
a = np.array(a)
b = [2*10**4,0,0]# m/s
b = np.array(b)
   

def minimum_initial_velocity(x_0,a,b,n):
    f_a = travel_to_the_moon(x_0,a)
    f_b = travel_to_the_moon(x_0,b)
    
    if (np.linalg.norm(a-b) < 0.1 and f_a*f_b < 0):
        return (a+b)/2
    if ( n > 100):
        print("Over 100 bisection iterations!")
        return -99

    if (travel_to_the_moon(x_0,(a+b)/2)*f_a < 0) :
        b = (a+b)/2
        n = n+1
        return(minimum_initial_velocity(x_0,a,b,n))
    
    else:
        a = (a+b)/2
        n = n+1
        return(minimum_initial_velocity(x_0,a,b,n))
        
v_0_min = minimum_initial_velocity(x_0,a,b,1)
print("The minimum initial speed to reach the moon is about", round(np.linalg.norm(v_0_min)*10**-3,2), "km/s.")
'''

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
