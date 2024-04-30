get_ipython().run_line_magic('matplotlib', 'inline')
import math
import matplotlib.pyplot as plot
import random
import copy
import numpy as np
import scipy.special as sp
import scipy.optimize as opt
import numpy as np

##########################################################################
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
    # print(t)S
    return a,b

###########################################################################

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

################################################################

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


################################################################


def fixpt(guess, err, func):#this function takes the guess and the permissible error as the input and finds the root using the fixed point method.
    x = guess
    x1 = x+3 # decalaring another variable for storage
    iter = 0 #to keep count of the number of iterations
    while abs(x - x1) >= err:
        x1 = x
        x = func(x)
        iter += 1
    return x, iter

###############################################################

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

def gaussian_quadrature(f, a, b, degree):
    nodes, weights =  np.polynomial.legendre.leggauss(degree)
    result = 0.5 * (b - a) * sum(w * f(0.5 * (b - a) * x + 0.5 * (a + b)) for x, w in zip(nodes, weights))
    return result

######################################################################

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

######################################################################

def secant_method(f, x0, x1, tol):
    t = 0
    h0 = f(x0)
    h1 = f(x1)
    
    while abs(h1) >= tol:
        x_new = x1 - h1 * (x1 - x0) / (h1 - h0)
        x0, x1 = x1, x_new
        h0, h1 = h1, f(x_new)
        t += 1

    #print("Root: ", x_new)
    #print("No. of iterations:", t)
    return x_new

######################################################################

def forward_euler(ode, y0, t0, t_end, dt):#ode is the derivative of y wrt t, y0, x0 are the initital boundary conditions, t_end is the end time, dt is the time skip
    t = np.arange(t0, t_end + dt, dt)
    y = np.zeros_like(t)
    y[0] = y0
    for i in range(1, len(t)):
        y[i] = y[i - 1] + dt * ode(t[i - 1], y[i - 1])

    return t, y



######################################################################

def backward_euler(ode, y0, t0, t_end, dt):
    t = np.arange(t0, t_end + dt, dt)
    y = np.zeros_like(t)
    y[0] = y0  

    for i in range (0, len(y)):
        def implicit_eq(y_new):
            return y_new - y[i - 1] - dt * ode(t[i], y_new)
        y[i] = opt.fsolve(implicit_eq, y[i - 1])[0]

    return t, y

######################################################################

def predictor_corrector(ode, y0, t0, t_end, dt):
    t = np.arange(t0, t_end + dt, dt)
    y = np.zeros_like(t)
    y[0] = y0

    for i in range(1, len(y)):
        y[i] = y[i-1] + dt*ode(t[i-1], y[i-1])
        y[i] = y[i-1]+dt*(ode(t[i-1], y[i-1])+ode(t[i], y[i]))/2
    return t, y



######################################################################


def RK_4(x, y, h, xf, func):# to  implement RK4 method to find the solution of a differential equation.
    x_val = [x]
    y_val = [y]
    x1 = x
    y1 = y
    while x < xf:
        k1 = h*func(x, y)
        k2 = h*func((x+0.5*h), (y+0.5*k1))
        k3 = h*func((x+0.5*h), (y+0.5*k2))
        k4 = h*func((x+h), y+k3)
        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x += (h)
        x_val.append(x)
        y_val.append(y)
    
    return x_val, y_val



######################################################################

def RK_4_multi(x_init, y_init, h, x_final, func):#to implement RK4 to solve n coupled ode. xinit, yinit, h, x_final are all first rank arrays of length n, function must take the array of x values and return the arrray of y values
    n = len(y_init)
    y_init = np.array(y_init)
    x_val = np.arange(x_init, x_final, h)
    m = len(x_val)
    y_val = np.zeros((n, m))
    x_val[0] = x_init
    y_val[:, 0] = y_init[:]
    for i in range (1, m):
        k1 = k2 = k3 = k4 = np.zeros(n)
        k1 = h*func(x_val[i-1], y_val[:, i-1])
        k2 = h*func((x_val[i-1]+0.5*h), (y_val[:, i-1]+0.5*k1))
        k3 = h*func((x_val[i-1]+0.5*h), (y_val[:, i-1]+0.5*k2))
        k4 = h*func((x_val[i-1]+h), y_val[:, i-1]+k3)
        y_val[:, i] = y_val[:, i-1]+(k1 + 2 * k2 + 2 * k3 + k4) / 6
    # print(x_val, y_val)
    
    return x_val, y_val



#################################################################################


def shooting (x, bv, func, guess1, guess2, tolerance = 0.001):#to implement shooting method to solve second order differential equation. x ocntains the initial and final values of x, bv are the boundary values at the initial and final values of x, funciton returns the differential equations. guess1 and guess2 are the two initila guess values.
    sol1 = RK_4_multi(x[0], [bv[0], guess1], h = 0.001, x_final=x[1], func = func)[1][0]
    sol2 = RK_4_multi(x[0], [bv[0], guess2], h = 0.001, x_final=x[1], func = func)[1][0]
    count = 0
    # print(sol2[-1]-bv[1])
    while abs(sol2[-1] - bv[1]) > tolerance:
        # print(sol2[-1])
        guess3 = guess2 + (guess1 - guess2)*(bv[1] - sol2[-1])/(sol1[-1] - sol2[-1])
        guess1 = guess2
        guess2 = guess3
        sol1 = RK_4_multi(x[0], [bv[0], guess1], h = 0.001, x_final=x[1], func = func)[1][0]
        sol2 = RK_4_multi(x[0], [bv[0], guess2], h = 0.001, x_final=x[1], func = func)[1][0]
        count+=1
    sol = RK_4_multi(x[0], [bv[0], guess2], h = 0.001, x_final=x[1], func = func)
    # print(count)
    return sol



#################################################################################




def leapfrog (pos_init, vel_init, acceleration, t0, t_fin, dt): #to implement the leapfrog algorithm to simulatte the dynamics of an object under force. Args - initial position, initial momentum, acceleration, initial and final time, and step size.
    #function returns the position array, velocuty arrayn and time stamps where the positions and velocitites are calculated...
    position = pos_init
    velocity = vel_init
    pos_arr = []
    vel_arr = [vel_init]
    velocity += acceleration(position)*(dt/2)
    time_arr_pos = np.arange(t0, t_fin, dt)
    time_arr_vel = np.arange(t0+(dt/2), t_fin-(dt/2), dt)

    for i in range (len(time_arr_vel)-1):
        vel_arr.append(velocity)
        pos_arr.append(position)
        position += velocity*dt
        velocity += acceleration(position)*(dt)
    pos_arr.append(position)
    vel_arr.append(velocity)
    position += velocity*dt
    velocity += acceleration(position)*(dt/2)
    pos_arr.append(position)
    vel_arr.append(velocity)
    time_arr_vel = np.append(t0, time_arr_vel)
    time_arr_vel = np.append(time_arr_vel, t_fin)
    return pos_arr, vel_arr, time_arr_pos, time_arr_vel




#################################################################################
#velocitry verlet

def verlet_method(f, y0, v0, t_span, h):
    num_steps = int((t_span[1] - t_span[0]) / h) + 1
    t_values = np.linspace(t_span[0], t_span[1], num_steps)
    
    y_values = np.zeros((num_steps, len(y0)))
    y_values[0] = y0
    
    # Initial step using Euler method
    y_temp = y0 + h * f(t_values[0], y0)
    y_values[1] = y0 + h * 0.5 * (f(t_values[0], y0) + f(t_values[1], y_temp))
    
    # Verlet method iterations
    for i in range(2, num_steps):
        y_values[i] = 2 * y_values[i-1] - y_values[i-2] + h**2 * f(t_values[i-1], y_values[i-1])
    
    return t_values, y_values

#VELOCITY VERLET

def velocity_verlet(f, x0, v0, t_span, h):
    
    num_steps = int((t_span[1] - t_span[0]) / h) + 1
    t_values = np.linspace(t_span[0], t_span[1], num_steps)
    x_values = np.zeros((num_steps, len(x0)))
    v_values = np.zeros((num_steps, len(v0)))

    x_values[0] = x0
    v_values[0] = v0

    for i in range(1, num_steps):
        # Verlet algorithm
        x_values[i] = x_values[i - 1] + h * v_values[i - 1] + 0.5 * h**2 * f(x_values[i - 1], t_values[i - 1])
        v_values[i] = v_values[i - 1] + 0.5 * h * (f(x_values[i], t_values[i - 1]) + f(x_values[i - 1], t_values[i - 1]))

    return t_values, x_values, v_values


#shooting



#boundary value 




#################################################################################

def finite_difference_bvp(f, a, b, alpha, beta, N):

    h = (b - a) / N
    x_values = np.linspace(a, b, N + 1)

    # Build the coefficient matrix and the right-hand side vector
    A = np.zeros((N - 1, N - 1))
    rhs = np.zeros(N - 1)

    for i in range(1, N):
        x_i = a + i * h
        A[i - 1, i - 1] = 2 + h**2 * f(x_i)
        rhs[i - 1] = -h**2 * f(x_i) * x_i

        if i < N - 1:
            A[i - 1, i] = -1
            A[i, i - 1] = -1

    # Adjust the system for boundary conditions
    rhs[0] -= alpha
    rhs[-1] -= beta

    # Solve the linear system
    y_interior = solve(A, rhs)

    # Combine interior and boundary values
    y_values = np.concatenate(([alpha], y_interior, [beta]))

    return x_values, y_values

#################################################################################




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


#################################################################################

    
def calculate_norm(r):
    """
    Calculates the Euclidean norm (L2 norm) of a vector r and returns the Euclidean norm of the vector r.
    """
    norm_squared = sum(element ** 2 for element in r)
    norm = norm_squared ** 0.5
    return norm


def conjugate_gradient(A, b, x0, tol=10**(-4), max_iter=100):
   
    x = np.array(x0, dtype=float)  # Initial guess
    r = b - np.dot(A, x)  # Initial residual
    p = r.copy()  # Initial search direction
    r_r = np.dot(r, r)
   
    
    for k in range(max_iter):
        Ap = np.dot(A, p)
        alpha = r_r / np.dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        
        r_r_next = np.dot(r, r)
        beta = r_r_next / r_r
        p = r + beta * p
        r_r = r_r_next
        
        if calculate_norm(r) < tol:
            break
        
    return x
  
    
#################################################################################

def initialize_grid(nx, ny):
    return np.zeros((nx, ny))

def poisson_equation(u, x, y, nx, ny, dx, dy, num_iterations):
    for _ in range(num_iterations):
        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                u[i, j] = 0.25 * (u[i+1, j] + u[i-1, j] + u[i, j+1] + u[i, j-1] - dx * dy * x[i, j] * np.exp(y[i, j]))

                
def set_boundary_conditions(u, x, y):
    # Set boundary conditions
    u[:, 0] = x[:, 0]  # u(x, 0) = x
    u[:, -1] = x[:, -1] * np.exp(1)  # u(x, 1) = xe
    u[0, :] = 0  # u(0, y) = 0
    u[-1, :] = 2 * np.exp(y[-1, :])  # u(2, y) = 2e^y

def save_output_table(u, x, y, nx, ny, filename):
    with open(filename, 'w') as file:
        file.write("x\t\ty\t\tu(x,y)\n")
        for i in range(nx):
            for j in range(ny):
                file.write(f"{x[i, j]:.2f}\t{y[i, j]:.2f}\t{u[i, j]:.6f}\n")


def plot_3d(x, y, u):
   
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(x, y)
    surf = ax.plot_surface(X, Y, u, cmap='magma', edgecolor='k', linewidth=0.5)
    
    # Add color gradient legend
    mappable = cm.ScalarMappable(cmap='magma')
    mappable.set_array(u)
    cbar = plt.colorbar(mappable, ax=ax)
    cbar.set_label('u(x, y)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u(x, y)')
    plt.title('3D Plot of the Solution to Poisson\'s Equation')
    plt.show()
    
################################################################################


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

def finite_difference_bvp(f, a, b, alpha, beta, N):

    h = (b - a) / N
    x_values = np.linspace(a, b, N + 1)
    print(x_values)
    # Build the coefficient matrix and the right-hand side vector
    A = np.zeros((N - 1, N - 1))
    rhs = np.zeros(N - 1)

    for i in range(1, N):
        x_i = a + i * h
        A[i - 1, i - 1] = 2 + h**2 * f(x_i)
        rhs[i - 1] = -h**2 * f(x_i) * x_i

        if i < N - 1:
            A[i - 1, i] = -1
            A[i, i - 1] = -1

    # Adjust the system for boundary conditions
    rhs[0] -= alpha
    rhs[-1] -= beta
    print(A)
    print(rhs)
    # Solve the linear system
    y_interior = seidel(A, rhs)
    print()
    # Combine interior and boundary values
    y_values = y_interior.insert(0,alpha)
    print(y_values)
    y_values = y_values.append(beta)
    return x_values, y_values


#####################################################################################

def crank_nicolson(M, alpha, T = 0.5, L = 8, k = 1):
    
    N = M**2
    x0, xL = 0, L
    dx = (xL - x0)/(M-1)

    # Final time 
    t0, tF = 0, T 

    # Time step
    dt = (tF - t0)/(N-1)

    # Stability parameter
    a0 = 1 + 2*alpha
    c0 = 1 - 2*alpha

    xspan = np.linspace(x0, xL, M)
    tspan = np.linspace(t0, tF, N)

    maindiag_a0 = a0*np.ones((1,M))
    offdiag_a0 = (-alpha)*np.ones((1, M-1))

    maindiag_c0 = c0*np.ones((1,M))
    offdiag_c0 = alpha*np.ones((1, M-1))

    #Left-hand side tri-diagonal matrix
    a = maindiag_a0.shape[1]
    diagonalsA = [maindiag_a0, offdiag_a0, offdiag_a0]
    A = sparse.diags(diagonalsA, [0,-1,1], shape=(a,a)).toarray()
    
    A[0,1] = (-2)*alpha
    A[M-1,M-2] = (-2)*alpha

    #Right-hand side tri-diagonal matrix
    c = maindiag_c0.shape[1]
    diagonalsC = [maindiag_c0, offdiag_c0, offdiag_c0]
    
    Arhs = sparse.diags(diagonalsC, [0,-1,1], shape=(c,c)).toarray()
    Arhs[0,1] = 2*alpha
    Arhs[M-1,M-2] = 2*alpha

    # ----- Initializes matrix U -----
    U = np.zeros((M, N))

    # Define the initial condition
    def u_init(x):
        return 4*x - x**2 / 2

    # Initial condition
    U[:,0] = u_init(xspan)
    
    f = np.arange(1, N+1) # Left Boundary Condition
    U[0,:] = 0
    f = U[0,:]
    
    g = np.arange(1, N+1) # Right Boundary Condition
    U[-1,:] = 0
    g = U[-1,:]
    
    for k in range(1, N):
        ins = np.zeros((M-2,1)).ravel()
        b1 = np.asarray([4*alpha*dx*f[k], 4*alpha*dx*g[k]])
        b1 = np.insert(b1, 1, ins)
        b2 = np.matmul(Arhs, np.array(U[0:M, k-1]))
        b = b1 + b2  # Right hand side
        U[0:M, k] = np.linalg.solve(A,b)  # Solve x=A\b 
                                          # for solving linear systems 
    return (U, tspan, xspan)





################################################################################

def My_Random(x0):
    a = 1103515245
    c = 12345
    m = 32768
    y = ((((a*(x0)+c))%m)/m)  
    return y



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


################################################################################

def Power_method_EV(A, num_iterations=1000, tol=1e-6):
    n = len(A)
    x = np.random.rand(n)  # Random initial guess for the eigenvector

    for _ in range(num_iterations):
        x1 = np.dot(A, x)
        eigenvalue = np.linalg.norm(x1)
        x1 = x1 / eigenvalue  # Normalize the eigenvector estimate
        if np.linalg.norm(x - x1) < tol:
            break
        x = x1

    return eigenvalue, x

################################################################################

def G_S_EV(A):
    n = len(A)
    Q = np.zeros((n, n))
    R = np.zeros((n, n))

    for j in range(n):
        v = A[:, j]
        for i in range(j):
            R[i, j] = np.dot(Q[:, i], A[:, j])
            v = v - R[i, j] * Q[:, i]
        R[j, j] = np.linalg.norm(v)
        Q[:, j] = v / R[j, j]

    return Q, R

def eignevalue(A, iterations = 50000):
    A1 = np.copy(A)
    # print(A_k)
    N = len(A)
    QQ = np.eye(N)
    for k in range(iterations):
        Q, R = G_S_EV(A1)
        A1 = R @ Q
        ev = []
    for i in range (0, A1.shape[0]):
        ev.append(A1[i, i])
    ev = np.array(ev)
    return ev, A1, QQ


################################################################################


def pol_fit(x_val, y_val, basis):
    avg = 0
    #print(x_val.shape)
    for x in (y_val):
        avg += x
    avg = avg/(x_val.shape)
    lhs = basis(x_val) @ basis(x_val).T
    rhs = basis(x_val) @ y_val.T
    par = np.linalg.inv(lhs)@rhs
    return par, np.linalg.cond(lhs)
################################################################################



def MLCG(a, m, c=0, x = 10):
	while True:
		x = (a * x + c) % m
		yield x / m      
################################################################################


def plot_MLCG(a,m,c,x):
    RN = []
    N = []
    for i in range(0,100):
        y = (a*x + c) % m
        RN.append(y/m)
        x = y
        N.append(i)  
    plot=plt.scatter(N, RN, label=f'a={a}, m={m}') 
    return plot
################################################################################



# Monte Carlo integration
def Monte_carlo(f, a, b, n, gen):
	sum = 0.0
	for _ in range(n):
		x = a + (b - a) * next(gen)  
        # Scale the random number to the interval [a, b]
		sum += f(x)
	return (b - a) * sum / n

################################################################################


# Monte Carlo integration with importance sampling
def monte_carlo_imp_sampling(f, p, inverse_cdf_p, n):
    samples = inverse_cdf_p(np.random.uniform(0, 1, n))
    weights = f(samples) / p(samples)
    return np.mean(weights), np.var(weights)




def product(a,b):
    a1 = list(a)
    b1 = list(b)
    g=[0,0,0,0,0]
    for i in range(0,5):
        sum1 = 0
        for j in range(0,5):
            sum1 = sum1+a1[i][j]*b1[j]
        g[i] = sum1
    return g    
            
def dot_product(a,b):
    a1 = list(a)
    b1 = list(b)
    sum1 = 0
    for i in range(0,5):
        sum1 = sum1 + a1[i]*b1[i]
    return sum1
            
            
           
            
            
#eigenvalue:

def eigen_value(A):
    x0 = [3,2,3,2,3]
    y = [3,2,3,2,3]
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
        # print(Ax0,"for",count,"iteration")
    sum2 = 0
    for i in range(0,len(Ax0)):
        sum2 = sum2 + Ax0[i]**2
    for i in range(0,len(Ax0)):
        Ax0[i] = Ax0[i]/(sum2**0.5)

    return Ax0,lambda1,count
    

