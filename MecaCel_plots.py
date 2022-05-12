import time
start_time = time.time()
import numpy as np
import matplotlib.pyplot as plt
import MecaCel_main as m


G = 2.959e-4                        
masses = [1, 0.000954786104043, 0.0002857214681]

tmax = 800 * 365 #jours
k = 90 #jours
n_k = int(tmax/k)

def resultat(methode, corps, masses):
    """
    methode : 
    0 pour Euler avant, 1 pour Heun, 2 pour RK4, 3 pour SV

    """
    if corps == 2:
        masses = masses[:2]
        q_Heun, p_Heun = m.start2(n_k,masses)
        if methode == 0:
            finalq, finalp = m.eula(q_Heun, p_Heun, k, n_k,masses)
        elif methode == 1:
            finalq, finalp = m.Heun(q_Heun, p_Heun, k, n_k,masses)
        elif methode == 2:
            finalq, finalp = m.RK4_new(q_Heun, p_Heun, k, n_k,masses)
        elif methode == 3:
            finalq, finalp = m.SV(q_Heun, p_Heun, k, n_k,masses)
        else:
            print("Les méthodes ne vont que de 0 à 3")
            finalq, finalp = 0, 0
    elif corps == 3:
        q_Heun, p_Heun = m.start3(n_k,masses)
        if methode == 0:
            finalq, finalp = m.eula(q_Heun, p_Heun, k, n_k,masses)
        elif methode == 1:
            finalq, finalp = m.Heun(q_Heun, p_Heun, k, n_k,masses)
        elif methode == 2:
            finalq, finalp = m.RK4_new(q_Heun, p_Heun, k, n_k,masses)
        elif methode == 3:
            finalq, finalp = m.SV(q_Heun, p_Heun, k, n_k,masses)
        else:
            print("Les méthodes ne vont que de 0 à 3")
            finalq, finalp = 0, 0
    else:
        print("Ce code ne gère que les cas à 2 et 3 corps")
        finalq, finalp = 0, 0
    return finalq, finalp, methode, masses
#------------------------------Simulation---------------------------------------
finalq, finalp, methode, masses_new = resultat(0, 2, masses)
#Premier argument => méthode : euler avant = 0, Heun = 1, RK4 = 2 et SV = 3
#Deuxieme argument => corps : 2 corps min. et 3 corps max.
#Troisieme argument => liste des masses des objets : onp eut laisser "masses" dans tous les cas de figures
#------------------------------Energie---------------------------------------

en = m.Energy(finalq, finalp, k, n_k,masses_new)
plt.plot(range(n_k), en)
plt.xlabel('Temps (jours)', fontsize = 15)
plt.ylabel('Énergie ($M_{s}*ua^{2}*jours^{-2}$)',fontsize = 15)
plt.title('Énergie',fontsize = 20)
if methode == 3:
    plt.xlim(0,250)
plt.show() 

#------------------------------Moment angulaire---------------------------------------

"""
plt.plot(range(n_k), m.angular_m(finalq, finalp, 0, n_k,masses_new),color = "Blue", label ='x')
plt.plot(range(n_k), m.angular_m(finalq, finalp, 1, n_k,masses_new),color = "Red", label ='y')
plt.plot(range(n_k), m.angular_m(finalq, finalp, 2, n_k,masses_new),color = "Green", label ='z')
plt.xlabel('Temps (jours)',fontsize = 15)
plt.ylabel('Moment angulaire ($M_{s}*ua^{2}*jours^{-1}$)',fontsize = 15)
plt.title('Moment angulaire SV',fontsize = 15)
plt.legend(loc="upper left")

plt.show()
"""
#------------------------------Orbite---------------------------------------
"""
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')

if len(masses_new) == 2 :
    ax.plot3D(finalq[:,3] - m.cm(finalq,0,n_k,masses_new),finalq[:,4] - m.cm(finalq,1,n_k,masses_new),finalq[:,5] - m.cm(finalq,2,n_k,masses_new), color = "Red", label ='Jupiter') #Jupiter
    ax.plot3D(finalq[:,0] - m.cm(finalq,0,n_k,masses_new),finalq[:,1] - m.cm(finalq,1,n_k,masses_new),finalq[:,2] - m.cm(finalq,2,n_k,masses_new), color = "Yellow", label = 'Soleil') #Soleil
if len(masses_new) == 3:
    ax.plot3D(finalq[:,6] - m.cm(finalq,0,n_k,masses_new),finalq[:,7] - m.cm(finalq,1,n_k,masses_new),finalq[:,8] - m.cm(finalq,2,n_k,masses_new), color = "Blue", label = 'Saturne') #Saturne
    ax.plot3D(finalq[:,3] - m.cm(finalq,0,n_k,masses_new),finalq[:,4] - m.cm(finalq,1,n_k,masses_new),finalq[:,5] - m.cm(finalq,2,n_k,masses_new), color = "Red", label ='Jupiter') #Jupiter
    ax.plot3D(finalq[:,0] - m.cm(finalq,0,n_k,masses_new),finalq[:,1] - m.cm(finalq,1,n_k,masses_new),finalq[:,2] - m.cm(finalq,2,n_k,masses_new), color = "Yellow", label = 'Soleil') #Soleil
ax.set_xlabel('x',fontsize = 15)
ax.set_ylabel('y',fontsize = 15)
ax.set_zlabel('z',fontsize = 15)
plt.legend(loc='upper left')
plt.title('Orbite (ua)')
plt.show()
print("--- %s seconds ---" % (time.time() - start_time))
"""