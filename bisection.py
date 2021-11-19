import numpy as np
   
def bisection(a,b,f,n):
    if (n > 10):
        print("Over 10 iterations!")
        return np.array([a,b])
    
    if ( f((a+b)/2)*f(a) < 0 ):
        b = (a+b)/2 
        n = n+1
        return(bisection(a,b,f,n))
    
    else:
        a = (a+b)/2
        n = n+1
        return(bisection(a,b,f,n))
    


