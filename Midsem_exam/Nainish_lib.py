import math
import numpy as np
import matplotlib.pylab as plt
from scipy import sparse



def LU_decomposition(z):
    LTriangle = []
    for y in range(len(z)):
        a = []
        for i in range(len(z)): #creating 2 zero matrix of equalm size
            a.append(0)
        LTriangle.append(a)

    UTriangle = []
    for y in range(len(z)):
        a = []
        for i in range(len(z)):
            a.append(0)
        UTriangle.append(a)


    for i in range(len(z)):                 
        for k in range(i, len(z)):    #decomposing matrix into L and U
            sum1 = 0                  # such that L*U = z
            for j in range(i):
                sum1 = sum1 + (LTriangle[i][j] * UTriangle[j][k])            
            UTriangle[i][k] = z[i][k] - sum1    
 
        for k in range(i, len(z)):
            if (i == k):
                LTriangle[i][i] = 1         #Doolittle LU to make diagonal terms 1 .                                
            else: 
                sum1 = 0                                                       
                for j in range(i):
                    sum1 = sum1 + (LTriangle[k][j] * UTriangle[j][i])
 
                LTriangle[k][i] = ((z[k][i] - sum1)/UTriangle[i][i])
    return LTriangle,UTriangle #returning L and U


############################################################

#Forward substitution:
def forward(z,b):
    i = 0
    X1 = []
    for k in range(0,len(z)): #creating a zero array
        X1.append(0)
    while i<len(z):
        j = 0
        ctr = 0
        while j<i:
            ctr = ctr + z[i][j]*X1[j]
            j = j+1
        X1[i] = (b[i]-ctr)/z[i][i]
        i = i+1
    return X1

############################################################

#Backward substitution:     
def backward(z,X1):
    i = len(z)-1
    X2 = []
    for l in range(len(z)):  #creating a zero array
        X2.append(0)
    while i>= 0:
        j = i
        sum1 = 0
        while j<len(z):
            sum1 = sum1 + z[i][j]*X2[j]
            j = j+1
            
        X2[i] = (X1[i]-sum1)/z[i][i]
        i = i - 1
    return X2

############################################################
def newton_raphson(a,z,f1,f1_dash):
    z=2.7
    e = 0.0000001 #Defining Tolerence
    d = 0.0000001
    k = 0
    while e < abs((a-z)) or abs(f1(a)) > d:
        z = a
        a = a - f1(a)/f1_dash(a)
        
        k = k+1
       
        print(a,"is the value for",k,'th iteration')
    return a


############################################################

def regulafalsi(a,b,f1):
    e = 0.000001 #Defining Tolerence
    d = 0.000001
    k = 0
    c = 9938927.369983709
    z = 0
    while e < abs((c-z)) or abs(f1(a)) > d:
        if f1(a)*f1(b) < 0:
            c = b - ((b-a)*f1(b))/(f1(b)-f1(a))
            z = c
            if f1(a)*f1(c) < 0:
                b = c
            else:
                a = c
        else:
            if abs(f1(a))>abs(f1(b)):
                a = a - 1.5*(b - a)
            else:
                b = b + 1.5*(b - a)
        print(a," is the root for" ,k ," iterations")
        k = k + 1 


############################################################
def bisection(a,b,f):
    if (f(a)*f(b)<=0):
        c = a
        while ((b-a)>= 0.000001):
            for i in range(15):
                c=(a+b)/2
                if (f(c)==0):
                    break
                if (f(c)*f(a)<0):
                    b=c
                else:
                    a=c
            
        return c
############################################################   
def bracket(a,b,f):
    t = 0
    while (f(a)*f(b)>=0):
        if abs(f(a))<abs(f(b)):
            
            b = b + 0.1
        else:
            a = a - 0.1
        bisection(a,b,f)
        t = t+1
    # print(t)
    return a,b

############################################################  
def simpson( a, b, N,f ):
    h = ( b - a )/N
    x = []
    fx = []
    i = 0
    while i<= N:
        x.append(a + i * h)
        fx.append(f(x[i]))
        i += 1
    sum1 = 0
    i = 0
    while i<= N:
        if i == 0 or i == N:
            sum1 = sum1 + fx[i]
        elif i % 2 != 0:
            sum1 = sum1 + 4 * fx[i]
        else:
            sum1 = sum1 + 2 * fx[i]
        i+= 1
    sum1 = round((sum1 * (h / 3)),9)
    return sum1
###################################################################





###################################################################

def rk4_li(t0, x0, v0, h, upper, x1,dvdt,dxdt):
    x = [x0]
    v = [v0]
    t = [t0]
    xi = x0
    vi = v0
    ti = t0
    zeta1 = v0
    while t0 <= upper:
        kx1 = h*(dxdt(vi,t0))
        kv1 = h*(dvdt(vi,xi,t0))
        
        kx2 = h*(dxdt(vi + (kv1/2), t0 + (h/2)))
        kv2 = h*(dvdt(xi + (kx1/2),vi + kv1/2, t0 +(h/2)))
        
        kx3 = h*(dxdt(vi+(kv2/2), t0+(h/2)))
        kv3 = h*(dvdt(xi+(kx2/2),vi+kv2/2,t0+(h/2)))
        
        kx4 = h*(dxdt((vi + kv3), (t0 + h)))
        kv4 = h*(dvdt((xi + kx3),vi+kv3,(t0 + h)))
             
        xi += (kx1 + 2*kx2 + 2*kx3 + kx4)/6
        vi += (kv1 + 2*kv2 + 2*kv3 + kv4)/6
        t0 += h
        
        
        x.append(xi)
        v.append(vi)
        t.append(t0)
    
    return x,v,t 

###################################################################
def interpolation_l(vl,vh,X1,V1,T1,X2,V2,T2,x1):
    n = len(X1) - 1
    if X1[len(X1)-1] < x1:
        if X2[len(X2)-1] > x1:
            dv = vh - vl
            k = vl + (dv/(X2[n]-X1[n]))*(x1 - X1[n])
    return k



def pde(l_x,l_t,N_x,N_t):
    import matplotlib.pylab as plt
    h_x = (l_x/N_x)
    h_t = (l_t/N_t)
    V0 = []
    X = []
    for i in range(0,N_x+1):
        # print(i)
        if h_x*i == 1:
            V0.append(300)
        else:
            V0.append(0)
        X.append(i)
    # print(V0)
    # print(0)
    X1 = []
    for i in range(0,len(X)):
        X1.append(X[i]/10)
    
    plt.plot(X1, V0)
    alpha  = (h_t/(h_x**2))
    # print(alpha)
    V1 = np.zeros(N_x+1)
    
    for j in range(0,1000):
        for i in range(0,N_x+1):
            if i == 0:
                V1[i] = (1-2*alpha)*V0[i] + alpha*V0[i+1]
            elif i == N_x:
                V1[i] = alpha*V0[i-1] + (1-2*alpha)*V0[i]
            else:
                V1[i] = alpha*V0[i-1] + (1-2*alpha)*V0[i] + alpha*V0[i+1]
        V0 = list(V1)
        
        if j==1 or j==2 or j==5 or j==10 or j==20 or j==50 or j==100 or j==200 or j==500 or j==1000:
            plt.plot(X1, V0)
    plt.xlabel("Length")
    plt.ylabel("Temperature")
    plt.title("Solution to 1d heat equation")
    plt.show()