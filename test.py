from numpy import *
from matplotlib.pyplot import *

def forward(k, Ts, T0, tend, dt):
    '''
    Solve T' = -k(T-Ts) when T(0) = T0 and t goes from 0 to t_end
    '''
    Nt = int(round(tend/float(dt)))
    tend = Nt*dt
    T = zeros(Nt+1)
    t = linspace(0,tend,Nt+1)
    T[0] = T0
    for n in range(0,Nt):
        T[n+1] = T[n]*(1-k*dt)+k*dt*Ts
    return t, T
dt = 0.5
k=1
Ts=273
T0 = 293
tend=10

t,T = forward(k, Ts, T0, tend,dt)

def exact(t,Ts, T0, k):
    return Ts + exp(-k*t)*(T0 - Ts)
    
def backward(k, Ts,T0, tend, dt):
    Nt = int(round(tend/float(dt)))
    T = zeros(Nt+1)
    tend = dt*Nt
    #t = linspace(0,tend,Nt+1)
    T[0] = T0
    for n in range(0,Nt):
        T[n+1] = (k*dt*Ts + T[n])/float(1+k*dt)
    return T
    
Tb = backward(k, Ts, T0, tend, dt)

def midpoint(k, Ts, T0, dt, tend):
    Nt = int(round(tend/float(dt)))
    tend = Nt*dt
    T = zeros(Nt+1)
    T[0] = T0
    for n in range(0,Nt):
        T[n+1] = (-k*dt*0.5*T[n]+T[n]+k*dt*Ts)/float(1+k*0.5*dt)
    return T

Tm = midpoint(k, Ts, T0, dt, tend)

t_e = linspace(0,tend,1001)
T_e = exact(t_e, Ts, T0, k)

plot(t_e, T_e, 'b-',
     t, T, 'r--o',
     t, Tb,'g--o',
     t, Tm,'y--o')
legend(['exact', 'forward', 'backward', 'midpoint'])
ylabel('Temperature [K]')
xlabel('time [s]')
title('Cooling')     
show()

