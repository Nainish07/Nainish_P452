import numpy as np
import math
import matplotlib.pyplot as plt

# MLCG
def MLCG(a, m, c=0, x = 10):
	while True:
		x = (a * x + c) % m
		yield x / m      

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


# Monte Carlo integration
def Monte_carlo(f, a, b, n, gen):
	sum = 0.0
	for _ in range(n):
		x = a + (b - a) * next(gen)  
        # Scale the random number to the interval [a, b]
		sum += f(x)
	return (b - a) * sum / n


# Monte Carlo integration with importance sampling
def monte_carlo_imp_sampling(f, p, inverse_cdf_p, n):
    samples = inverse_cdf_p(np.random.uniform(0, 1, n))
    weights = f(samples) / p(samples)
    return np.mean(weights), np.var(weights)
