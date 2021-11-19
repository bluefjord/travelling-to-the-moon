import numpy as np
import pandas as pd
from random import  uniform

G = 6.7 * 10 ** (-11)
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
    def move(self,dt):
        a = self.acceleration(self.pos)
        self.pos = self.pos + self.v *dt + 0.5*a*dt**2
        self.v = self.v +a*dt

def reachTheMoon(myrocket,moon,earth,dt):
    t = 0

    while t < 1*10**6:

        if myrocket.v < 0:
            return False
        
        myrocket.move(dt)

        t = t+dt
        
        if myrocket.pos > moon.pos - moon.radius:
            return True

    return False

def timeTaken(myrocket,moon,earth,dt):
    t = 0

    while t < 1*10**6:

        if myrocket.v <0:
            return -99
  
        myrocket.move(dt)

        t = t+dt
        
        if myrocket.pos > moon.pos - moon.radius:
            return t



def randomise(variable, maxError):
    variable = uniform(variable+maxError, variable-maxError)
    return variable

