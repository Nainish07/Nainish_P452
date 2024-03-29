import math
import numpy as np
import matplotlib.pyplot as plt




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

###########################################################
def prnt_matrix(a):
  for i in range(0,len(a)):
    for j in range(0,len(a[i])):
      print(a[i][j],end = " ")
    print("")
  print("")

##########################################################





        
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
        print('for x',(i+1),end =" root is equal to ") # for printing the last element of each row.
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
        
    return z


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
  
    
    
#QUESTION 4 ASSIGNMENT 2

def matrix_A_ij(x):
    m = 0.2
    N = len(x)
    delta = 1.0 / N
    result = np.zeros_like(x)
    
    for i in range(N):
        result[i] += (delta + m) * x[i]
        result[i] -= 2 * delta * x[i]
        result[i] += delta * x[(i + 1) % N]  # Periodic boundary condition
        result[i] += delta * x[(i - 1) % N]  # Periodic boundary condition
        result[i] += m ** 2 * delta * x[i] 
    #print(result)
    return result

def conjugate_fly(matrix_A_ij, b, x0, tol=10**(-6), max_iter=100):
    it = 0
    x = x0.copy()  # Initial guess
    r = b - matrix_A_ij(x)  # Initial residual
    p = r.copy()  # Initial search direction
    residue_norms = [calculate_norm(r)]  # List to store residue norms

    for k in range(max_iter):
        Ap = matrix_A_ij(p)
        alpha = np.dot(r, r) / np.dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        
        beta = np.dot(r, r) / np.dot(r - alpha * Ap, r - alpha * Ap)
        p = r - alpha * Ap + beta * p

        residue_norm = calculate_norm(r)
        residue_norms.append(residue_norm)
        it= it+1
        if residue_norm < tol:
            break

    return x, residue_norms, it


def conjugate_inv(matrix_A_ij, b, x0, tol=10**(-6), max_iter=100):
    N = len(b)
    inverse_columns = []
    
    for i in range(N):
        # Create the right-hand side vector for solving Ax = e_i
        ei = np.zeros(N)
        ei[i] = 1
        
        # Solve the equation Ax = e_i using Conjugate Gradient method
        x, _, _ = conjugate_fly(matrix_A_ij, ei, x0, tol, max_iter)
        
        # Append the solution (column of the inverse matrix) to the list
        inverse_columns.append(x)
    
    # Stack the columns of the inverse matrix horizontally to form the complete inverse matrix
    A_inv = np.column_stack(np.round(inverse_columns,4))
    return A_inv