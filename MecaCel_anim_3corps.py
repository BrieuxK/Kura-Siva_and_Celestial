import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import MecaCel_main as m



tmax = 2000 * 365 #jours
k = 90 #jours
n_k = int(tmax/k)
masses = [1, 0.000954786104043, 0.0002857214681]

q_Heun, p_Heun = m.start3(n_k,masses)

finalq, finalp = m.SV(q_Heun, p_Heun, k, n_k,masses)

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1, projection='3d')

x= np.concatenate([finalq[:,3], finalq[:,6]])
y= np.concatenate([finalq[:,4], finalq[:,7]])
z= np.concatenate([finalq[:,5], finalq[:,8]])

points, = ax.plot(x, y, z, 'o', c = 'r', markersize = 3) #Ca doit être Jupiter celle là
dots, = ax.plot(x, y, z, 'o', c = 'b', markersize = 3)  #Peut être Saturne
line1, = ax.plot(x, y, z, '-', c = 'b') #Ligne pour Saturne
line2, = ax.plot(x, y, z, '-', c = 'r') #Ligne pour Jupiter
txt = fig.suptitle('')
plt.grid(False)

def update_points(num, x, y, z, points, dots, line1, line2):
    txt.set_text("nbr d'années={:d}".format(num*k//365))

    new_x = x[num] ####################"
    new_y = y[num]
    new_z = z[num]
    
    new_x3 = x[:num]
    new_y3 = y[:num]
    new_z3 = z[:num]
    
    new_x1 = x[num+n_k]
    new_y1 = y[num+n_k]
    new_z1 = z[num+n_k]
    
    new_x2 = x[n_k:num+n_k]
    new_y2 = y[n_k:num+n_k]
    new_z2 = z[n_k:num+n_k]

    points.set_data(new_x,new_y)
    points.set_3d_properties(new_z, 'z')
    dots.set_data(new_x1,new_y1)
    dots.set_3d_properties(new_z1, 'z')
    line1.set_data(new_x2,new_y2)
    line1.set_3d_properties(new_z2, 'z')
    line2.set_data(new_x3,new_y3)
    line2.set_3d_properties(new_z3, 'z')

    return points,dots,line1,line2,txt

ax.set_xlim3d([-10.0, 10.0])

ax.set_ylim3d([-10.0, 10.0])

ax.set_zlim3d([-1, 1])


ax.plot(0, 0, marker = 'o', markersize=15, color="yellow")
ani=animation.FuncAnimation(fig, update_points, frames=n_k//3, fargs=(x, y, z, points, dots, line1, line2))
writergif = animation.PillowWriter(fps=15)
ani.save('3corpsSV.gif',writer=writergif,progress_callback=lambda i, n: print(f'Saving frame {i} of {n}'))
plt.show()