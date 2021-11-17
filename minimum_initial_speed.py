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
       # print("a",a)
    # print(self.pos)
        self.pos = self.pos + self.v *dt + 0.5*a*dt**2
     #   print(self.v)
        self.v = self.v +a*dt
     #   print(self.v)

def reachTheMoon(myrocket):
    t = 0

    while t < 1*10**6:

       # print(myrocket.pos,myrocket.v)

        if myrocket.v < 0:
            #print("Didn't reach the moon")
            return False
        
        myrocket.move()

        t = t+dt
        
        if myrocket.pos > moon.pos - moon.radius:
     #       print("true")
            return True

    return False

'''
    def move(self):
        
        a_n = self.acceleration(self.pos)
        a_n1 = self.acceleration(self.pos + dr)

        self.v[0] = (self.v[0] + np.sqrt(self.v[0]**2 + 4*( (a_n ) )*dr))/2
        self.v[1] = (self.v[1] + np.sqrt(self.v[1]**2 + 4*( (a_n1) )*dr))/2
        self.pos = self.pos + dr
'''

'''
    def getTimeIncrease(self):

        dt_small = (-3*10**9)*np.ones(2)
        dt_large = (3*10**9)*np.ones(2)

        for i in range(0,2):

            if self.pos < zeroAccelerationPoint and self.pos+dr < zeroAccelerationPoint:

                v_n = self.v[i]
                #a_n = self.a
                a_n1 = self.acceleration(self.pos+dr)
 
                dt_small[i] = dr / v_n
                dt_large[i] = ( -v_n + np.sqrt(v_n**2 + 4 * a_n1 * dr ) ) / ( 2 * a_n1 ) 
                #print(dt_small, dt_large)

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
'''


'''
def timeTaken(x_0,v_0):

    bodies = [earth,moon]
    myrocket = rocket(x_0,v_0, bodies)

    t_smaller = 0
    t_larger = 0

    counter = 0

    while counter < 1*10**10:
        
        distancerocketearth = abs(myrocket.pos - earth.pos)
        distancerocketmoon = abs(myrocket.pos - moon.pos)

       # print(myrocket.v)
        if distancerocketearth > 2 * moon.pos:
            break
            #return -1
        if distancerocketearth < earth.radius:
            break
            #return -2
        if distancerocketmoon < moon.radius:
            break


        else:
            dt = myrocket.getTimeIncrease()
            if np.isnan(dt[0]):
                print("going back to earth")
                break
            
            t_smaller = t_smaller + dt[0]
            t_larger = t_larger + dt[1]
            myrocket.move()
            counter = counter +1


    return [t_smaller, t_larger, myrocket,counter]


'''


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
