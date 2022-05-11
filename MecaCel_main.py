import numpy as np


G = 2.959e-4                        #1.989*1e30   #Masse Soleil
#masses = [1, 0.000954786104043, 0.0002857214681]

def start2(n_k,masses): #Pour 2 corps
    q_Heun = np.zeros((n_k, 2 * 3))  #1st col : XSun, 2nd : YSun, 3rd : ZSun,...., 6th : ZJup
    p_Heun = np.zeros((n_k, 2 * 3))
    q_Heun[0] = np.array([0, 0, 0,4.6580839119399,-1.7945574704081,-0.0967627432056])
    p_Heun[0] = np.array([0, 0, 0,0.0026255239604*masses[1],0.0074041879836*masses[1],-0.0000894925649*masses[1]])
    return q_Heun, p_Heun

def start3(n_k,masses): #Pour 3 corps
    q_Heun = np.zeros((n_k, len(masses) * 3))  #1st col : XSun, 2nd : YSun, 3rd : ZSun,...., 6th : ZJup
    p_Heun = np.zeros((n_k, len(masses) * 3))

    #COnditions initiales dans la 1ère ligne de chaque matrice
    q_Heun[0] = np.array([0, 0, 0,4.6580839119399,-1.7945574704081,-0.0967627432056,6.9600794880566,-7.0666123077577,-0.1541266614081])
    p_Heun[0] = np.array([0, 0, 0,0.0026255239604*masses[1],0.0074041879836*masses[1],-0.0000894925649*masses[1],0.0036673089657*masses[2],0.0039101598352*masses[2],-0.0002138877692*masses[2]])
    
    return q_Heun, p_Heun

def dq(pos, imp, masses): 
    zq = np.zeros(3*len(masses))
    for i in range(len(masses)):
        j = i*3
        zq[j] = imp[j]/masses[i]
        zq[j+1] = imp[j+1]/masses[i]
        zq[j+2] = imp[j+2]/masses[i]
    return zq


def dp(param_pos,param_imp,masses):
    zp = np.zeros(len(param_pos))
    for i in range(len(masses)):
        j = i*3
        for m in range(len(masses)):
            n = m*3
            if i != m :
                zp[j]  += -2*((G*masses[i]*masses[m])/distance([param_pos[j],param_pos[j+1],param_pos[j+2]], [param_pos[n],param_pos[n+1],param_pos[n+2]])**3) * (param_pos[j]-param_pos[n])
                zp[j+1]+= -2*((G*masses[i]*masses[m])/distance([param_pos[j],param_pos[j+1],param_pos[j+2]], [param_pos[n],param_pos[n+1],param_pos[n+2]])**3) * (param_pos[j+1]-param_pos[n+1])
                zp[j+2]+= -2*((G*masses[i]*masses[m])/distance([param_pos[j],param_pos[j+1],param_pos[j+2]], [param_pos[n],param_pos[n+1],param_pos[n+2]])**3) * (param_pos[j+2]-param_pos[n+2])
    return zp*0.5

def distance(pos1, pos2):
    """
    pos1,pos2 = [X,Y,Z]
    """
    dist = np.sqrt(((pos1[0]-pos2[0])**2)+((pos1[1]-pos2[1])**2)+((pos1[2]-pos2[2])**2))
    return dist

def F_Heun(position, impulsion, k,masses):
    """
    position,impulsion : array [x1 y1 z1 x2 y2 z2]
    """
    q1 = k * dq(position,impulsion,masses)
    #print("q1 :",q1)
    p1 = k * dp(position,impulsion,masses)
    q_tilde = position + q1
    #print("q_tilde :", q_tilde)
    p_tilde = impulsion + p1
    q2 = k * dq(q_tilde, p_tilde,masses)
    p2 = k * dp(q_tilde, p_tilde,masses)
    return 0.5*(q1+q2), 0.5*(p1+p2)

def Heun(q_Heun, p_Heun, k, n_k,masses):
    """
    q/p_Heun : matrice (n_k,6)
    CHaque colonne représente une donnée 
    Chaque ligne représente un pas de temps
    """
    for i in range(n_k-1):
        delta_q, delta_p = F_Heun(q_Heun[i],p_Heun[i], k,masses)
        q_Heun[i+1] = q_Heun[i] + delta_q
        p_Heun[i+1] = p_Heun[i] + delta_p
    return q_Heun, p_Heun

def eula(pos,imp, k, n_k,masses):
    for i in range(n_k-1):
        pos[i+1] = pos[i] + k * dq(pos[i],imp[i],masses)
        imp[i+1] = imp[i] + k * dp(pos[i],imp[i],masses)
    return pos, imp

def RK4_new(q, p, k, n_k,masses):
    for i in range(n_k-1):  
        k1q, k1p = k * dq(q[i], p[i],masses), k * dp(q[i], p[i],masses)
        k2q, k2p = k * dq(q[i]+k1q/2., p[i]+k1p/2.,masses), k * dp(q[i]+k1q/2., p[i]+k1p/2.,masses)
        k3q, k3p = k * dq(q[i]+k2q/2., p[i]+k2p/2.,masses), k * dp(q[i]+k2q/2., p[i]+k2p/2.,masses)
        k4q, k4p = k * dq(q[i]+k3q, p[i]+k3p,masses), k * dp(q[i]+k3q, p[i]+k3p,masses)  #k2, k3?
        ###
        q[i+1] = q[i] + (k1q+2*k2q+2*k3q+k4q)/6
        p[i+1] = p[i] + (k1p+2*k2p+2*k3p+k4p)/6
    return q,p

def SV(q_SV,p_SV, k, n_k,masses):
    for i in range(n_k-1):
        p_mid = p_SV[i] + 0.5 * k * dp(q_SV[i], p_SV[i],masses)
        q_SV[i+1] = q_SV[i] + k * dq(q_SV[i], p_mid,masses)
        p_SV[i+1] = p_mid + 0.5 * k * dp(q_SV[i+1],p_mid,masses)
    return q_SV,p_SV

def Energy(result_q, result_p, k, n_k,masses): 
    energies = np.zeros(n_k)
    if len(masses) == 3:
        for i in range(n_k):       
            T_Sun = (result_p[i,0]**2 + result_p[i,1]**2 + result_p[i,2]**2)/(2*masses[0])
            T_Jup = (result_p[i,3]**2 + result_p[i,4]**2 + result_p[i,5]**2)/(2*masses[1])
            T_Sat = (result_p[i,6]**2 + result_p[i,7]**2 + result_p[i,8]**2)/(2*masses[2])
            V_SoJ = G*masses[0]*masses[1]/distance(result_q[i,:3],result_q[i,3:6])
            V_SoSat = G*masses[0]*masses[2]/distance(result_q[i,:3],result_q[i,6:])
            V_JSat = G*masses[1]*masses[2]/distance(result_q[i,3:6],result_q[i,6:])
            energies[i] = (T_Jup + T_Sun + T_Sat) - V_SoJ - V_SoSat - V_JSat
    if len(masses) == 2:
        for i in range(n_k):       
            T_Sun = (result_p[i,0]**2 + result_p[i,1]**2 + result_p[i,2]**2)/2*masses[0]
            T_Jup = (result_p[i,3]**2 + result_p[i,4]**2 + result_p[i,5]**2)/2*masses[1]
            V = G*masses[0]*masses[1]/distance(result_q[i,:3],result_q[i,3:])
            #print(T_Jup)
            energies[i] = (T_Jup + T_Sun) - V
    return energies

def angular_m(result_q, result_p, coord, n_k,masses): #Moment angulaire
    """
    result_q/_p => (n_k,6) : matrice de solutions finales
    coord => 0 pour x, 1 pour y, 2 pour z
    """
    ang_m = np.zeros(n_k)
    if len(masses) == 3:
        ang_m = result_q[:,0 + coord] * result_p[:,0 + coord] + result_q[:,3 + coord] * result_p[:,3 + coord] + result_q[:,6 + coord] * result_p[:,6 + coord]
    if len(masses) == 2:
        ang_m = result_q[:,0 + coord] * result_p[:,0 + coord] + result_q[:,3 + coord] * result_p[:,3 + coord]
    """
    return : moment angulaire par rapport à une seule coordonnée
    """
    return ang_m

def cm(result_q, coord, n_k,masses):  #Centre de masse
    """
    coord => 0 pour x, 1 pour y, 2 pour z
    """
    cmass = np.zeros(n_k)
    if len(masses) == 3:
        cmass = 1 * result_q[:,0 + coord] + 0.000954786104043 * result_q[:,3 + coord] + 0.0002857214681 * result_q[:,6 + coord]
        m_tot = 1 + 0.000954786104043 + 0.0002857214681
    if len(masses) == 2:
        cmass = 1 * result_q[:,0 + coord] + 0.000954786104043 * result_q[:,3 + coord]
        m_tot = 1 + 0.000954786104043
    """
    return : centre de masse par rapport à une seule coordonnée
    """
    return cmass/m_tot
