import numpy as np
import pandas as pd
import matplotlib as plt
from random import  uniform


class planet:

    def __init__(self, pos, radius, mass):
        self.pos = pos  # m
        self.radius = radius  # m
        self.mass = mass  # kg

G = 6.7 * 10 ** (-11)

p_e = 0
r_e = 6 * 10 ** 6
m_e = 6 * 10 ** 24

p_m = 4 * 10 ** 8
r_m = 2 * 10 ** 6
m_m = 1 * 10 ** 23

dr = 1*10**4

h = 100*10**3

zeroAccelerationPoint = ( p_m ) / ( np.sqrt(m_e / m_m) + 1 )

earth = planet(p_e, r_e, m_e)
moon = planet(p_m, r_m, m_m)


class rocket:

    def __init__(self, pos, v, bodies):
    
        self.pos = pos
        self.v = v
        self.bodies = bodies

    def acceleration(self,pos):
        a = 0
        for body in self.bodies:
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

            #else:
            #    print(self.pos < zeroAccelerationPoint and self.pos+dr < zeroAccelerationPoint)
            #    return [-99,-99]

        return [min(dt_small), max(dt_large)]
            
            
'''
    def t_1(self):
        if self.v != 0:
            t_1 = dr/self.v
        else:
            t_1 = -99
        return t_1

    def t_2(self,pastZeroAccerlationPoint):
        
        self.a = self.acceleration()
        
        
        
        discriminant = self.v**2 + 4 * self.a * dr 

        if discriminant < 0:
            print("Discrimant less than zero, initial velocity too slow")
            return -99

        t_2a = ( -self.v + np.sqrt(discriminant) ) / ( 2 * self.a ) 
        t_2b = ( -self.v - np.sqrt(discriminant) ) / ( 2 * self.a ) 
            
        if pastZeroAccerlationPoint: 
           t_2 = t_2a 

        else: 
           t_2 = max(t_2a, t_2b)
           
        return t_2

'''

def timeTaken(x_0,v_0):

    bodies = [earth,moon]
    myRocket = rocket(x_0,v_0, bodies)

    t_smaller = 0
    t_larger = 0

    counter = 0

    while counter < 1*10**10:
        
        distanceRocketEarth = abs(myRocket.pos - earth.pos)
        distanceRocketMoon = abs(myRocket.pos - moon.pos)

        if distanceRocketEarth > 2 * moon.pos:
            break
            #return -1
        if distanceRocketEarth < earth.radius:
            break
            #return -2
        if distanceRocketMoon < moon.radius:
            break


        else:
            dT = myRocket.getTimeIncrease()
            if dT[0] < 0 or dT[1] < 0:
                print("Going backwards in time")
                print(dT)
                break
            
            t_smaller = t_smaller + dT[0]
            t_larger = t_larger + dT[1]
            myRocket.move()
            counter = counter +1


    return [t_smaller, t_larger, myRocket,counter]
'''
        else:
            myRocket.move()
            
            if myRocket.pos <= zeroAccelerationPoint:
                if myRocket.t_1() == -99:
                    break
                if myRocket.t_2(False) == -99:
                    break
                t_smaller = t_smaller + myRocket.t_1()
                t_larger = t_larger + myRocket.t_2(False)

                
            else:
                #if counter % 1000== 0:
                #    print(myRocket.a)
                t_larger = t_larger + myRocket.t_1()
                t_smaller = t_smaller + myRocket.t_2(True)

            counter = counter + 1
'''

T = timeTaken(h + r_e,np.array([6.8*10**4, 6.8*10**4]))

print("The zero acceleration point is", zeroAccelerationPoint)
print("The final position of the rocket", T[2].pos)
print("The counter value:", T[3])

print("The smallest estimate for the time take to reach the mooon: ",T[0] , " and the largest:", T[1])
print("The difference between the smallest and the largest is:", T[1] - T[0])

'''
def timeTaken(v_0):
    t = 0
    bodies = [earth,moon]
    myRocket = rocket(x_0,v_0,bodies)

    while t < 1*10**8:
        
        distanceRocketEarth = np.linalg.norm(myRocket.pos - earth.pos)
        distanceRocketMoon = np.linalg.norm(myRocket.pos - moon.pos)

        if distanceRocketEarth > 2 * apogee:
            return -1
        if distanceRocketEarth < earth.radius:
            return -2
        if distanceRocketMoon < moon.radius:
            return -3

        else:
            my_rocket.move()
            t = t + dt
'''

a = 1000*10**3 # m/s
b = 1300*10**3 # m/s


def randomise(variable, maxError):
    variable = uniform(variable+maxError, variable-maxError)
    return variable

N=10

v_0_array = np.linspace(12*10**3, 15*10**3, num=N)

time_mean_array = np.zeros(N)
time_sd_array   = np.zeros(N)

for i in range(N):

    time = np.zeros((N,2))
    for j in range(N):
        p_m = randomise(0.4055 * 10 ** 9, 0.00005 * 10 ** 9)
        p_e = 0
    
        m_m = randomise(0.07346 * 10 ** 24, 0.000005 * 10 * 24)
        m_e = randomise(5.9724 * 10 ** 24, 0.00005 * 10 ** 24)
    
        r_m = randomise( 1736.0 * 10 ** 3, 0.05*10**3)
        r_e = randomise(6356.752 * 10 ** 3, 0.0005 * 10 ** 3)
    
    
        
        zeroAccelerationPoint = ( p_m ) / ( np.sqrt(m_e / m_m) + 1 )
    
        earth = planet(p_e, r_e, m_e)
        moon = planet(p_m, r_m, m_m)
    
    
        
        T = timeTaken(h + r_e, np.array([v_0_array[i], v_0_array[i]]))
        time[j] = [T[0], T[1]]
    
    #print("Here are the times for an initial velocity of", v_0, ":", time)
    
    
    '''
    # Let's calculate the standard deviation of the time and the mean time

    var = 0
    mean = 0
    n = 10
    for i in range(n):
        var += ( ( (time[i][1] - time[i][0]) / 2)**2 )
        mean += (time[i][1] +  time[i][0])/2
    
    mean = mean/n
    var = var/n
    '''


    time_mean = np.mean(time)
    time_mean_array[i] = time_mean

    #time_max = np.max(time)
    #time_min=np.min(time)

    #time_error = 0.5*(time_max-time_min)
    time_sd = np.std(time)
    time_sd_array[i] = time_sd

    
    #print("The error bounds of the time data for", v_0," m/s is between ", time_max, "and",time_min," seconds")

'''
    print("The mean of the time data for", v_0_array[i], " m/s is", time_mean," seconds")
    print("So the result is for ", v_0_array[i], "m/s is", time_mean, "+-", time_sd,"seconds")


'''




d = {'vel': v_0_array, 'mean time': time_mean_array, 'std time': time_sd_array}

#time_vel_data= pd.DataFrame(v_0_array,time_mean_array, time_sd_array)

time_vel_data= pd.DataFrame(d)
time_vel_data.to_csv("time_vel_data.csv")



