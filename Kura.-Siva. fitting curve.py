#Courbe passant par les différentes valeurs de longueur critiques, déterminées pour des viscosités différentes.

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt;

N = 1024
dt = 0.05
T = 200
n_t = int(T/dt)

t = np.linspace(0,T,n_t)
k = np.arange(-N/2, N/2)

def KS(L, visc, N = 1024, T = 200):
   
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

L_all = []
Lcrit_list  = []

def lcrit(visc):
   
    Lcrit_list  = []
   
    for i in range(visc):
        print(i)
        A_all = []
        for L in range(2,50):
           
            x, t, u, u_spec2 = KS(L,i)

            A2 = sum(u[-1]**2)*1/len(u)*1/L
            A = np.sqrt(A2)
            A_all.append(A)
       
        for l in range(len(A_all)):
       
            if A_all[l]<0.001:
                continue
            else:
                Lcrit_list.append(l+2)
                break
               
    return Lcrit_list[1:]

res = lcrit(12)
print(res)

t = np.arange(1, 12)
plt.plot(t, res, linestyle="None", marker="o")
plt.plot(t, 7*t**((np.log(13/7))/(np.log(4))), color = "Black", label ="approx.")
#plt.savefig("longcrit plot")
#plt.show()

r = np.polyfit(np.sqrt(t), res, 1) #[5.56843937 1.40166647]
print(r)

y = []
for i in t:
    ordo = r[0] * np.sqrt(i) + r[1]
    y.append(ordo)

plt.plot(t,y, color ='Red', label ="fit")
plt.xlabel("Viscosité",fontsize = 10)
plt.ylabel("Longueur critique",fontsize = 10)
plt.legend()
plt.show()
