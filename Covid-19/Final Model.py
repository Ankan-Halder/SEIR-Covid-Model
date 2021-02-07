#SEIR without vital dynamics
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt




N  = input("Please enter the total population: ")
I0 = input("Please enter the initial number of infected: ")
R0 = input("Please enter the initial number of removed or recovered: ")
E0 = input("Please enter the initial number of exposed: ")

S0 = N - E0 - I0 - R0 




# Contact rate:  beta
# recovery rate: gamma
# exposed rate:  sigma
# direct recovery rate: meu

beta  = input("Please enter the contact rate: ")
gamma = input("Please enter the recovery rate : ") 
sigma = input("Please enter the exposed rate : ") 
meu   = input("Please enter the direct recovery rate : ") 

n = input("Please enter the  number of days: ")

t = np.linspace(0, n, num = n)




# The SIR model differential equations.

def SEIR_modal(y,t,N,beta, gamma,sigma,meu):
    S, E, I, R = y
    dSdt =  - beta * S * I/N 
    dEdt = beta * S * I/N  - sigma * E - meu * E
    dIdt = sigma * E  - gamma * I 
    dRdt = gamma * I + meu * E

    return dSdt, dEdt, dIdt, dRdt

y0 = S0, E0, I0, R0

sol = odeint(SEIR_modal, y0, t, args=(N,beta, gamma,sigma,meu))
S, E, I, R = sol.T




# Plot the data on four separate curves for S(t),E(t), I(t) and R(t)

fig = plt.figure(facecolor='w')
plt.plot(t, S/N, 'b',  label='Susceptible')
plt.plot(t, E/N, 'y',  label='Exposed')
plt.plot(t, I/N, 'r',  label='Infected')
plt.plot(t, R/N, 'g',  label='Recovered with immunity')
plt.title("SEIR Covid-Model")
plt.xlabel("Days")
plt.ylabel("Fraction of Population")
plt.ylim(0,.03)
plt.legend(loc="upper right")
plt.grid()
plt.show()