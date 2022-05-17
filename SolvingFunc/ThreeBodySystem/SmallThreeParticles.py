from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
from scipy.integrate import odeint
matplotlib.rcParams.update({'font.size': 18})

def gravity3body(y,t,G,m1,m2,m3):
    # momentum and position
    px1,px2,px3,py1,py2,py3,qx1,qx2,qx3,qy1,qy2,qy3 = y
    pp12 = -1*G*m1*m2/(((qx1-qx2)**2+(qy1-qy2)**2)**(3/2))
    pp13 = -1*G*m1*m3/(((qx1-qx3)**2+(qy1-qy3)**2)**(3/2))
    pp23 = -1*G*m2*m3/(((qx2-qx3)**2+(qy2-qy3)**2)**(3/2))
    def qq(p,m):
        return p/m
    def ppdiff(qa,qb):
        return qa-qb
    dydt = []
    dydt.append(pp12*ppdiff(qx1,qx2)+pp13*ppdiff(qx1,qx3)) # px1'
    dydt.append(pp12*ppdiff(qx2,qx1)+pp23*ppdiff(qx2,qx3)) # px2'
    dydt.append(pp13*ppdiff(qx3,qx1)+pp23*ppdiff(qx3,qx2)) # px3'
    dydt.append(pp12*ppdiff(qy1,qy2)+pp13*ppdiff(qy1,qy3)) # py1'
    dydt.append(pp12*ppdiff(qy2,qy1)+pp23*ppdiff(qy2,qy3)) # py2'
    dydt.append(pp13*ppdiff(qy3,qy1)+pp23*ppdiff(qy3,qy2)) # py3'
    dydt.append(qq(px1,m1)) # qx1'
    dydt.append(qq(px2,m2)) # qx2'
    dydt.append(qq(px3,m3)) # qx3'
    dydt.append(qq(py1,m1)) # qy1'
    dydt.append(qq(py2,m2)) # qy2'
    dydt.append(qq(py3,m3)) # qy3'
    return dydt

### start main script ###
G = 1
m1 = 2
m2 = 2
m3 = 2

### assume three particles have momentum (1,1), (0,1) and (0,-1)###
### assume three particles initially at (0,0), (1,0) and (-1,0)###
y0 = [1.,0.,0.,1.,1.,-1.,0.,1.,-1.,0.,0.,0]

### define the time grid ###
dt = 0.01
t = np.arange(0.0,10.0,dt)

### generate the solution ###
sol = odeint(gravity3body,y0,t,args=(G,m1,m2,m3))
posx1,posy1 = sol[:,6],sol[:,9]
posx2,posy2 = sol[:,7],sol[:,10]
posx3,posy3 = sol[:,8],sol[:,11]

'''
### create the animation frames ###   
iframe = 0
for i in range(len(t)):
    ### output the current frame ###  
    output_filename = '%d.png' % iframe   
    fig = plt.figure(figsize=(6,6))
    plt.scatter(posx1[i],posy1[i],c='r',s=m1*65)
    plt.scatter(posx2[i],posy2[i],c='g',s=m2*65)
    plt.scatter(posx3[i],posy3[i],c='b',s=m3*65) 
    plt.tick_params(
        axis='both',       # changes apply to both axes
        which='both',      # both major and minor ticks are affected
        bottom=True,      # ticks along the bottom edge are off
        top=True,         # ticks along the top edge are off
        right=True,       # ticks along the right edge are off
        left=True,        # ticks along the left edge are off
        labelbottom=True, # labels along the bottom edge are off
        labelleft=True    # labels along the left edge are off
    )
    plt.grid()
    plt.axis([-2,6.5,-0.75,2.5])
    plt.title('Simulation Time = '+str(int(i*dt)))
    plt.savefig(output_filename,format="png")
    plt.close(fig)
    iframe = iframe+1

### stitch together the frames to make a movie ### 
os.system("rm hw05b.mp4")
os.system("ffmpeg -i %d.png -vf scale=432x432 hw05b.mp4")
os.system("rm *.png") 
'''

### plot the solution ###
fig,ax=plt.subplots(figsize=(6,6))
ax.plot(posx1,posy1,'r',label='m1=%.2f'%(m1))
ax.plot(posx2,posy2,'g',label='m2=%.2f'%(m2))
ax.plot(posx3,posy3,'b',label='m3=%.2f'%(m3))
ax.legend()
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.grid()
#plt.savefig('traj.png',format='png')
plt.show()

### plot the trajectory with time ###
fig = plt.figure(figsize=[6,6])
ax = plt.axes(projection='3d')
ax.scatter(posx1,posy1,t,s=1,label='m1=%.2f'%(m1))
ax.scatter(posx2,posy2,t,s=1,label='m2=%.2f'%(m2))
ax.scatter(posx3,posy3,t,s=1,label='m3=%.2f'%(m3))
ax.set(xlabel='x',ylabel='y',zlabel='t')
ax.legend()
#plt.savefig('traj3D.png',format='png')
plt.show()