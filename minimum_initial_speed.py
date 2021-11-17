import numpy as np
import pandas as pd
from random import  uniform



G = 6.7 * 10 ** (-11)
dr = 1*10**6
dt = 1
h = 100*10**3

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

    def acceleration(self,pos):
        a = 0
        for body in self.bodies:
            r = body.pos - self.pos
            r_abs = abs(r)
            a = a + ((G * body.mass) / (r_abs) ** 3) * r
            
        return a
    def move(self):
        a = self.acceleration(self.pos)
        self.pos = self.pos + self.v *dt + 0.5*a*dt**2
        self.v = self.v +a*dt

def reachTheMoon(myrocket):
    t = 0

    while t < 1*10**6:

        if myrocket.v < 0:
            return False
        
        myrocket.move()

        t = t+dt
        
        if myrocket.pos > moon.pos - moon.radius:
            return True

    return False

def randomise(variable, maxError):
    variable = uniform(variable+maxError, variable-maxError)
    return variable

bins = 10
rep = 10

v_0_set = np.linspace(10.9869*10**4, 10.9871*10**4, num=bins)

v_0_count = np.zeros(bins)


for i in range(bins):
    
    for j in range(rep):
        p_m = randomise(0.4055 * 10 ** 9, 0.00005 * 10 ** 9)
        p_e = 0
        
        m_m = randomise(0.07346 * 10 ** 24, 0.000005 * 10 * 24)
        m_e = randomise(5.9724 * 10 ** 24, 0.00005 * 10 ** 24)
        
        r_m = randomise( 1736.0 * 10 ** 3, 0.05*10**3)
        r_e = randomise(6356.752 * 10 ** 3, 0.0005 * 10 ** 3)


        earth = planet(p_e, r_e, m_e)
        moon = planet(p_m, r_m, m_m)

        bodies = [earth,moon]

        myrocket = rocket(h, v_0_set[i], bodies)
        #print(myrocket.pos,myrocket.v)
 
        if(reachTheMoon(myrocket)):
            v_0_count[i] = v_0_count[i] +1


d = {'bins': v_0_set, 'count': v_0_count}
v_0_distribution_data = pd.DataFrame(d)


v_0_distribution_data.to_csv("v_0_distribution_data.csv")

print(v_0_distribution_data)

#print("Some estimates for the minimum initial velocity", minimum_v_array)
