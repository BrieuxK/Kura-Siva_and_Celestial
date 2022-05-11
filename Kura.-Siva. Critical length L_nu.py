#Diagramme de la longueur critique en fonction de L pour une viscosité donnée.

import numpy as np
import matplotlib.pyplot as plt

N = 1024
visc = 2
dt = 0.05
T = 300
n_t = int(T/dt)

t = np.linspace(0,T,n_t)
k = np.arange(-N/2, N/2)

def KS(L, N = 1024, T = 200):
    x = np.linspace(0,L - L/N ,N)
    
    #Conditions initiales
    CI = np.cos(2*np.pi*x/L) + 0.1*np.cos(4*np.pi*x/L)
    
    #Matrices de solutions, réelle et esp.spectral
    u = np.zeros((n_t,N),dtype = complex) #40K lignes, 1024 colonnes
    u_spec1 = np.zeros((n_t, N),dtype = complex)
    u_spec2 = np.zeros((n_t, N),dtype = complex) #Chloé => F(u²) =/= [F(u)]² 
                                                 #Donc on crée une matrice pour F(u) et une autre pour F(u²)
    
    #On remplit les matrices avec les C.I.
    u[0] = CI
    u_spec1[0] = 1/N * np.fft.fftshift(np.fft.fft(CI))
    u_spec2[0] =  1/N * np.fft.fftshift(np.fft.fft(CI**2))
    
    #Transfo de Fourier de la partie lin.
    f_L = (2*np.pi*k/L)**2 - visc*(2*np.pi*k/L)**4
    
    #Pour la dérivée
    a = 1j * 2*np.pi/L * k
    
    #On remplit le reste des matrices
    u_spec1[1] = ((1 + dt*0.5*f_L)/(1 - dt*0.5*f_L))*u_spec1[0] -  a/2*(3/2*(u_spec2[0]) - 1/2*(u_spec2[0])) * dt/(1 - dt*0.5*f_L)
    u[1] = N * np.fft.ifft(np.fft.ifftshift(u_spec1[1]))
    u_spec2[1] = 1/N * np.fft.fftshift(np.fft.fft(u[1]**2))
    for i in range(1,n_t - 1):
        u_spec1[i+1] = (1 + dt*0.5*f_L)/(1 - dt*0.5*f_L)*u_spec1[i] - a/2*(3/2*(u_spec2[i]) - 1/2*(u_spec2[i-1])) * dt/(1 - dt*0.5*f_L)
        #fftshift
        u[i+1] = N * np.fft.ifft(np.fft.ifftshift(u_spec1[i+1]))
        u_spec2[i+1] = 1/N * np.fft.fftshift(np.fft.fft(u[i+1]**2))
    
    return x,t,u,u_spec2

A_all = []
L_all = []

L_crit = 0

for L in range(2,20):
    L_all.append(L)    
    x, t, u, u_spec2 = KS(L)
    print(L)
    """
    A2 = np.trapz(u**2,x)  #méthode des trapèzes => info BAC2 
    A = np.sqrt(A2)
    A_all.append(A)
    """
    
    A2 = sum(u[-1]**2)/len(u)*1/L
    A = np.sqrt(A2)
    A_all.append(A)
    
    if A > 10e-3 and L_crit == 0:
        L_crit = L
        print("valeur crit",L_crit)

#print(np.shape(A_all[2]),"test")

plt.plot(L_all,A_all)
plt.xlabel('L')
plt.ylabel('A')
plt.title(f"Critical length for visc = {visc}")
plt.show()
