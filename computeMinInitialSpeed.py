import bisection
import numpy as np
import voyageToTheMoon as voyage


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

#myrocket = voyage.rocket(h+r_e, v_0_set[i], bodies)

a = 10*10**3
b = 16*10**3

def f(v):
    myrocket = voyage.rocket(h+r_e, v, bodies)
    
    result = 2*voyage.reachTheMoon(myrocket,moon,earth,dt)-1

    return result
ans = bisection.bisection(a,b,f,1)

print("The minimum initial velocity is between (km/s)", ans/1000)
 

# Let's calculate the time taken at 10 percent
# the minimum initial speed

vTenPercent = 1.1*np.mean(ans)

myrocket = voyage.rocket(h+r_e, vTenPercent, bodies)
time = voyage.timeTaken(myrocket,moon,earth,dt)

print("The time taken to reach the moon at 10 percent the minimum initial speed",f'{vTenPercent/100:.3g}' ,"(km/s) is", f'{time/3600:.3g}', "hours")



