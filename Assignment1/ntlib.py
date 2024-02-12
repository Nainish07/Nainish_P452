import math
import numpy as np
import matplotlib.pylab as plt


# def My_Random(x0):
#     a = 1103515245
#     c = 12345
#     m = 32768
#     y = (((a*(x0)+c))%m)    
#     return y

def My_Random(x0):
    a = 1103515245
    c = 12345
    m = 32768
    y = ((((a*(x0)+c))%m)/m)  
    return y

############################################################





def matrix(n,m,x):
  z = []
  x_1 = []
  for k in range(0,len(x)):
   x_1.append(x.pop())

  for i in range(0,n):
    a = []
    for j in range(0,m):
      a.append(x_1.pop())
    z.append(a)
  return z

############################################################

# printin of a matrix
def prnt_matrix(a):
  for i in range(0,len(a)):
    for j in range(0,len(a[i])):
      print(a[i][j],end = " ")
    print("")
  print("")

############################################################

def Product_1(a,b):
  N=[]
  
  for i in range(0,len(a)):
    h = []
    for j in range(0,len(b[0])):
      s = 0
      for k in range(0,len(b)):
        temp = a[i][k]*b[k][j]
        s = temp + s
      h.append(round(s,3))
    N.append(h)
  prnt_matrix(N)
  print(" ")

        

############################################################


def Product_2(a,b): # dot product
  z = []
  s = 0
  for i in range(0,len(a)):
    s = s +  a[i][0]*b[i][0]
  z.append(s)
  print(z)
  print(" ")


############################################################

# printin of a matrix
def prnt_matrix(a):
  for i in range(0,len(a)):
    for j in range(0,len(a[i])):
      print(a[i][j],end = " ")
    print("")
  print("")

############################################################

def max_swap(z,a):
  for i in range(a,len(z)): #for swapping(sorting) the rows with largest element and smallest element 
    for j in range(i+1,len(z)):
      if abs(z[i][a]) < abs(z[j][a]):
        m = z[i]
        z[i] = z[j]
        z[j] = m
      else:
        break
############################################################  


def norm_row(z,a):
  x = []  #for creating All the diagonal elements 1
  for i in range(0,len(z)+1):
    c = z[a][i]
    x.append(c/z[a][a])
  z[a] = x
  
############################################################

def coloumn(z,a):
  for i in range(a+1,len(z)):
    x = []
    if abs(z[i][a]) > 0:  #for making a elements 0 that are below the diagonal elemnt
      for j in range(0,len(z[i])):
        c = z[i][j] - z[a][j]*z[i][a]
        x.append(c)
      z[i] = x

############################################################
  

def reverse_coloumn(z,a): #for making a elements 0 that are above the diagonal elemnt
  for i in range(a-1,-1,-1):
    if abs(z[i][a]) > 0:
      for j in range(len(z[i])-1,-1,-1):
        c = z[i][j]-z[a][j]*z[i][a]
        z[i][j] = c
  return z

############################################################

def give_solutions(z):
    for i in range(0,len(z)): # the output matrix is an augumented matrix 
        print((i+1),end =" root is equal to ") # for printing the last element of each row.
        print(z[i][len(z)])

############################################################
        
# for gauss-jordan elemination
def gj_elimination(z):
  for i in range(0,len(z)): 
    max_swap(z,i)  #for creating a lower triangular matrix 
    norm_row(z,i) 
    coloumn(z,i)
  for j in range(len(z)-1,-1,-1):
    reverse_coloumn(z,j) #for creating all other elements above diagonal 0.
  give_solutions(z)
  return z


############################################################

def create_augument(z,x):
    for i in range(0,len(z)):
        z[i].append(x[i])
        print(z[i])
    return z


############################################################

def seidel(a,b):
  L = []
  v = []
  for i in range(0,len(a)):
    L.append(0) #creating 2 zero row matrix that has our guess solution 0.
    v.append(0)
  # print(L)
  # print(v)
  k = 0
  e = 10
  while e > 0.000001: #setting the epsilon 
    for i in range(0,len(a)):
      sum1 = 0
      sum2 = 0
      for j in range(0,i): 
        sum1 = sum1 + a[i-1][j-1]*L[j-1] 
      for j in range(i+1,len(a)):
        sum2 = sum2 + a[i-1][j-1]*L[j-1]
      v[i-1] = L[i-1] # storing the values of previous iteration
      L[i-1] = (b[i-1] - sum1 - sum2)/a[i-1][i-1]
    e = (abs(L[3])-abs(v[3]))/abs(L[3])  # calculation of episilon after each iteration
    k = k + 1
    print(L)
    print(k)
    if k > 30:  #stopping the loop after 30 iterations
      print('diverging')
      break

############################################################
        
def jacobi(z,b):
    X =  []
    v = []
    for i in range(0,len(z)):
        X.append(0)
        v.append(0) # creating 2 zero arrays of same length
    # print(X)
    k = 0
    e = 10
    while e > 0.001 : #setting the episilon 
        for i in range(0,len(z)):
            v[i] = X[i] #storing the values of previous iteration
            sum1 = 0
            for j in range(0,len(z)):
                if j!=i:
                    sum1 = sum1 + z[i][j]*X[j]
            X[i] = (b[i] - sum1)/z[i][i] 
        e = (abs(X[4])-abs(v[4]))/abs(X[4]) #calculationg the epsilon
        k=k+1
        if k > 500: #reduce while using in exam
            print('diverging')
            break
    return X
        

############################################################       
        
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


def is_symetric(z):
    for i in range(0,len(z)):
        for j in range(0,len(z)):
            if z[i][j] != z[j][i]:
                return False
    return True

############################################################           
    
def Transpose_sq(z):
    x = []
    for i in range(0,len(z)):
        y = []
        for j in range(0,len(z)):
            y.append(z[j][i])
        x.append(y)
    return x
 
############################################################           

def cholesky(z):
    if is_symetric(z):                                       
        i=0  #checking symitric matrix
        while i <len(z):
            j=0
            sum1=0
            while j<i:
                sum1 = sum1 + z[j][i]*z[j][i]
                j = j + 1
            z[i][i]=(z[i][i]-sum1)**(0.5)                       
            j=i+1
            while j<len(z):
                k=0
                sum1=0
                while k<i:
                    sum1 = sum1 + z[i][k]*z[k][j]
                    k = k + 1
                z[j][i]=(z[j][i]-sum1)/z[i][i]                  
                z[i][j]=z[j][i]
                j=j+1
            i=i+1
        i=0
        while i <len(z):                                        
            j=i+1
            while j<len(z): #making all the elements above the diagonal 0
                z[i][j]=0   #To get a trangular matrix
                j=j+1
            i=i+1 
    return z


############################################################


def deflation(c,a):
    if len(c) != 1:
        c[1] = c[1] + c[0]*a
        for i in range(2,len(c)):
            c[i] = c[i] + a*c[i - 1]
        c.pop()
    else:
        print("cannot deflate")
    return c

############################################################


def copy_list(c):
    cn = []
    for i in range(0,len(c)):
        a = c[i]
        cn.append(a)
    return cn

############################################################


def derivative(z):
    
    for i in range(0,len(z)):
        z[i] = z[i]*(len(z)-i-1)
    z.pop()
    
    return z

############################################################    

def fn(c,z):
    sum1 =0
    for i in range(0,len(c)):
        sum1 = sum1 + c[i]*(z**(len(c)-i-1))
    return sum1

############################################################


def lagurre(c,z,degree):
    t = 0
    f = 7
     
    k1 = copy_list(c)
    k2 = copy_list(c)
    k2 = derivative(k2)
    
    while t < degree:
        if abs(fn(c,z))<0.0000001:
            print(fn(c,z))
            print(z, "is a root")
        else:
            k1 = derivative(k1)
            k2 = derivative(k2)
            count = 0
            while abs(fn(c,z))>0.0000001:
                f = z
                
                sum1 = 0
                for i in range(0,len(k1)):
                    sum1 = sum1 + k1[i]*(z**(len(k1)-i-1))
                G = sum1/fn(c,z)
               
                sum2 = 0
                for i in range(0,len(k2)):
                    sum2 = sum2 + k2[i]*(z**(len(k2)-1-i))
                H = G**2 - (sum2/fn(c,z)) 
                
                n = len(c)-1
                if G < 0:
                    a = n/(G - ((n-1)*(n*H - G**2))**0.5)
                else:
                    a = n/(G + ((n-1)*(n*H - G**2))**0.5)
                z = z - a
                count = count + 1
        print(z,"is a root")
        # print(count,'is the number of iterations used to get',z)
        c = deflation(c,z)
        k1 = copy_list(c)
        k2 = copy_list(c)
        k2 = derivative(k2)
        t=t+1
        # print(t,"round of while loop")
    # print((-1)*c[1],"is a root")

############################################################    

def interpolation(x,y,c):
    sum = 0
    for i in range(0,len(x)):
        product  = 1
        for k in range(0,len(x)):
            if k != i:
                product = product * ((c-x[k])/(x[i]-x[k])) 
        #print(product)
        d = product
        sum = sum + d*y[i]
    print(sum)

############################################################

def least_sq_fit(X,Y,degree):
    # plt.plot((X),Y)
    # plt.show()  
    V= []
    for i in range(0,degree):
        a = []
        for j in range(0,degree):
            a.append(0)
        V.append(a)
    z = []
    c = []
    for j in range(0,2*len(V)-1):
        sum1= 0
        sum2 = 0
        for i in range(0,len(X)):
            sum1 =sum1+(X[i]**j)
            sum2 = sum2 + Y[i]*(X[i]**(j))
        z.append(sum1)
        c.append(sum2)
    c.pop()
    # c.pop()
    for i in range(0,len(V)):
        for j in range(0,len(V)):
            V[i][j] = z[i+j]
    return V,c

############################################################

def newton_raphson(a,f1,f1_dash):
    z = 10
    e = 0.0000001
    d = 0.0000001
    k = 0
    while e < abs((a-z)) or abs(f1(a)) > d:
        z = a
        a = a - f1(a)/f1_dash(a)
        k = k+1
    return(a)
    print(a,"is the value for",k,'th iteration')
    # # print("")
    # # print(a,"is the root")

############################################################
def regulafalsi(a,b,f1):
    e = 0.000001
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

def find_root(a,b,f1):
    e = 0.000001
    d = 0.000001
    k = 0
    while e < abs((b-a)) or abs(f1(a)) > d:
        if f1(a)*f1(b) < 0:
            c = (a+b)/2
            if f1(c)*f1(a) < 0:
                b = c
            else:
                a = c
        else:
            if abs(f1(a))>abs(f1(b)):
                a = a - 1.5*(b - a)
            else:
                b = b + 1.5*(b - a)
        k = k + 1 
        
        print(a,"is the root for ",k,"iterations")
        


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
    
def bracket(a,b,f):
    t = 0
    while (f(a)*f(b)>=0):
        if abs(f(a))<abs(f(b)):
            a = a - 0.1
        else:
            b = b + 0.1
        bisection(a,b,f)
        t = t+1
    # print(t)
    return a,b




###################################################################


def diagonal_dominant(z,b):
    i=0
    n=len(z)
    for i in range(n):               
        m=z[0][i]              
        for j in range(n):
            if z[j][i]>m:
                m=z[j][i]
                z[i],z[j]=z[j],z[i]
                b[i],b[j]=b[j],b[i]
        
                if z[j][j]==0:
                    z[i],z[j]=z[j],z[i]
                    b[i],b[j]=b[j],b[i]
    return z,b


###################################################################

def linear_fit(X,Y,sigma):
    # plt.scatter((X),Y)
    # plt.show()  
    sum_x = 0
    sum_y = 0
    sum_x2 = 0
    sum_xy = 0
    sum_y2 = 0
    N = len(X)
    for i in range(0,len(X)):
        sum_x = sum_x + (X[i]/sigma[i])
        sum_x2 = sum_x2 + ((X[i]**2)/sigma[i])
        sum_y = sum_y + (Y[i]/sigma[i])
        sum_xy = sum_xy + ((X[i]*Y[i])/sigma[i])
        sum_y2 = sum_y2 + ((Y[i]**2)/sigma[i])
    a1 = (N*sum_xy-sum_x*sum_y)/(N*sum_x2-(sum_x**2))
    a2 = (sum_y - a1*sum_x)/(N)
    
    delta_x = N*sum_x2-(sum_x**2)
    delta_y = N*sum_y2-(sum_y**2)
    r = ((N*sum_xy - sum_x*sum_y)**2)/(delta_x*delta_y)
    print(a1,"is the slope")
    print(a2, "is the intercept")
    print(r,"is the preasons r^2 value")
    # return a1,a2,r


################################################################### 
    


def midpoint_method(a,b,N,f):
    h = (abs(a - b))/N
    x = []
    y = []
    # print(h)
    for i in range(1,N+1):
        x.append((2*a + (2*i-1)*h)/2)
        
        y.append(f(x[i-1]))
    # print(y)
    # print(x)
    sum1 = 0
    for j in range(0,len(y)):
        sum1 = sum1 + y[j]*h
    sum1 = round(sum1,9)
    return(sum1)

######################################################################

def trapezoidal(a,b,N,f):
    h = (abs(a - b))/N
    x = []
    y = []
    # print(h)
    
    for i in range(0,N+1):
        x.append((a + i*h))
        y.append(f(x[i]))
    sum1 = 0
    for i in range(1,len(x)):
        sum1 = sum1 + (h/2)*(y[i-1]+y[i])
    sum1 = round(sum1,9)
    return(sum1)


######################################################################


def monte_carlo(a,b,x_0,N,f):
    X = []
    Y = []
    sum_y2 = 0
    sum_y = 0
    sum_x = 0
    for i in range(0,N):
        x = My_Random(x_0)
        X.append(x*(b-a)+a)
        y = f(x)
        Y.append(y)
        x_0 = x
        sum_y = sum_y+y
        sum_y2 = sum_y2 + y**2
        
    sigma_sq = (sum_y2)/N -(sum_y/N)**2
    f_N = ((b-a)/N)*(sum_y)
    return(f_N)
    

###################################################################


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


#######################################################################################




def product(a,b):
    a1 = list(a)
    b1 = list(b)
    g=[0,0,0]
    for i in range(0,3):
        sum1 = 0
        for j in range(0,3):
            sum1 = sum1+a1[i][j]*b1[j]
        g[i] = sum1
    return g    
            
def dot_product(a,b):
    a1 = list(a)
    b1 = list(b)
    sum1 = 0
    for i in range(0,3):
        sum1 = sum1 + a1[i]*b1[i]
    return sum1
            
            
            
            
            
#eigenvalue:

def eigen_value(A):
    x0 = [3,2,3]
    y = [3,2,3]
    k = 6
    lambda1 = 0
    Ax0 = product(A,x0)
 
    Ax0x0 = dot_product(product(A,x0),x0)
    denominator = Ax0x0
    count = 0
    while abs(lambda1 - k) > 0.001:
        k = lambda1
        numerator = dot_product(product(A,Ax0),x0)
        lambda1 = (numerator/denominator)
        Ax0 = product(A,Ax0)
        denominator = numerator 
        count = count + 1
        print(Ax0,"for",count,"iteration")
    sum2 = 0
    for i in range(0,len(Ax0)):
        sum2 = sum2 + Ax0[i]**2
    for i in range(0,len(Ax0)):
        Ax0[i] = Ax0[i]/(sum2**0.5)

    print(Ax0)
    print(lambda1)
    print(count)
    
    
    
######################################################################################



def shm(t0,x0,v0,h,upper,dvdt,dxdt):
    x = [x0]
    v = [v0]
    t = [t0]
    xi = x0
    vi = v0
    ti = t0
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


#######################################################################################

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
    plt.plot(X, V0)
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
            plt.plot(X, V0)
    plt.show()
    
##########################################################################################

def interpolation_l(vl,vh,X1,V1,T1,X2,V2,T2,x1):
    n = len(X1) - 1
    if X1[len(X1)-1] < x1:
        if X2[len(X2)-1] > x1:
            dv = vh - vl
            k = vl + (dv/(X2[n]-X1[n]))*(x1 - X1[n])
    return k
            
    
###################################################################################




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

######################################################################################################


def pde2(l_x,l_t,N_x,N_t):
    h_x = (l_x/N_x)
    h_t = (l_t/N_t)
    V0 = []
    X = []
    for i in range(0,N_x+1):
        # print(i)
        V0.append(300*abs(math.sin(math.pi*h_x*i*0.5)))
        X.append(i)
    # print(V0)
    # print(0)
    plt.plot(X, V0)
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
        if j==10 or j==20 or j==50 or j==100 or j==200 or j==500 or j==1000:
            plt.plot(X, V0)
    plt.show()
    
#######################################################################################


def sho_2(d2ydx2, dydx, x0, y0, z0, xf, st):
   
    x = [x0]
    y = [y0]
    z = [z0]      # dy/dx

    n = int((xf-x0)/st)     # no. of steps
    for i in range(n):
        x.append(x[i] + st)
        k1 = st * dydx(x[i], y[i], z[i])
        l1 = st * d2ydx2(x[i], y[i], z[i])
        k2 = st * dydx(x[i] + st/2, y[i] + k1/2, z[i] + l1/2)
        l2 = st * d2ydx2(x[i] + st/2, y[i] + k1/2, z[i] + l1/2)
        k3 = st * dydx(x[i] + st/2, y[i] + k2/2, z[i] + l2/2)
        l3 = st * d2ydx2(x[i] + st/2, y[i] + k2/2, z[i] + l2/2)
        k4 = st * dydx(x[i] + st, y[i] + k3, z[i] + l3)
        l4 = st * d2ydx2(x[i] + st, y[i] + k3, z[i] + l3)

        y.append(y[i] + (k1 + 2*k2 + 2*k3 + k4)/6)
        z.append(z[i] + (l1 + 2*l2 + 2*l3 + l4)/6)

    return x, y, z

#####################################################################################################333