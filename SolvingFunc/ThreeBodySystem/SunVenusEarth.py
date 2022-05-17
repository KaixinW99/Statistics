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
G = 6.67e-11
m1 = 1988500e24
m2 = 4.87e24
m3 = 5.97e24

### assume three particles initially at (0,0), (-108.56e9,0) and (149.72e9,0)###
### assume three particles have momentum (0,0), (0,-m2*35e3) and (0,m3*29.8e3)###
y0 = [0.,0.,0.,0.,-m2*35e3,m3*29.8e3,0.,-108.56e9,149.72e9,0.,0.,0]

### define the time grid ###
### in 100 seconds range ###
dt = 86400
t = np.arange(0.0,365.25*86400.0,dt)

### generate the solution ###
sol = odeint(gravity3body,y0,t,args=(G,m1,m2,m3))
posx1,posy1 = sol[:,6],sol[:,9]
posx2,posy2 = sol[:,7],sol[:,10]
posx3,posy3 = sol[:,8],sol[:,11]

### create the animation frames ###   
iframe = 0
for i in range(len(t)):
    ### output the current frame ###  
    output_filename = '%d.png' % iframe   
    fig = plt.figure(figsize=(6,6))
    plt.scatter(posx1[i],posy1[i],c='r',s=432690/100)
    plt.scatter(posx2[i],posy2[i],c='g',s=3760.4/100)
    plt.scatter(posx3[i],posy3[i],c='b',s=3958.5/100) 
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
    plt.axis([-1.6e11,1.6e11,-1.6e11,1.6e11])
    plt.title('Simulation Time = '+str(int(i*dt)))
    plt.savefig(output_filename,format="png")
    plt.close(fig)
    iframe = iframe+1

### stitch together the frames to make a movie ### 
os.system("rm sun.mp4")
os.system("ffmpeg -i %d.png -vf scale=432x432 sun.mp4")
os.system("rm *.png") 

'''
### plot the solution ###
plt.figure(figsize=(10,10))
plt.plot(posx1,posy1,'b',label='m1=%.2f'%(m1))
plt.plot(posx2,posy2,'g',label='m2=%.2f'%(m2))
plt.plot(posx3,posy3,'r',label='m3=%.2f'%(m3))
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.show()
'''
