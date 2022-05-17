############## Choice 1 ##############
# import package
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
matplotlib.rcParams.update({'font.size': 16})

### animation parameters ###
xmin =  0.0
xmax = 10.0
ymin =  0.0
ymax = 10.0
tmin =  0.0
tmax = 60.0
vmin =  0.0
vmax =  1.0
dt   =  0.05

### define particle(s) ###
### (nb - ideally particles should move << 1 radius in dt) ### 
npart = 20
px = np.random.uniform(low=xmin+1.,high=xmax-1.,size=npart)
py = np.random.uniform(low=ymin+1.,high=ymax-1.,size=npart)
vx = np.random.uniform(low=vmin,high=vmax,size=npart)
vy = np.random.uniform(low=vmin,high=vmax,size=npart)

r = np.random.uniform(0.05,0.15,npart)  # particle radius in simulation units
m = 1.0*np.power(r,3)   # mass of one particle
rp = 65*r              # particle radius in points(?)

g = -1.0 # gavity constant
ax = np.full_like(px,0.0)
ay = np.full_like(px,g)

colors  = ['b']*npart

def collision(px1,px2,py1,py2,r1,r2,vx1,vx2,vy1,vy2,m1,m2):
    p1 = np.array([px1,py1],dtype=float)
    p2 = np.array([px2,py2],dtype=float)
    v1 = np.array([vx1,vy1],dtype=float)
    v2 = np.array([vx2,vy2],dtype=float)
    if (px1-px2)**2+(py1-py2)**2 <= (r1+r2)**2:
        p1 -= v1*dt   ### the collision takes 1 time interval (avoiding two balls sticking)###
        p2 -= v2*dt   ### the collision takes 1 time interval (avoiding two balls sticking)###
        const = 2*np.dot(v1-v2,p1-p2)/np.dot(p1-p2,p1-p2)/(m1+m2)
        v1 -= const*m2*(p1-p2)
        v2 -= const*m1*(p2-p1)    
    return v1[0],v1[1],v2[0],v2[1],p1[0],p1[1],p2[0],p2[1]

### create the animation frames ###   
iframe = 0

for t in np.arange(tmin,tmax,dt):
    ### output the current frame ###  
    output_filename = '%d.png' % iframe   
    fig = plt.figure(figsize=(6,6))
    plt.scatter(px,py,c=colors,s=rp*rp)   
    plt.tick_params(
        axis='both',       # changes apply to both axes
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        right=False,       # ticks along the right edge are off
        left=False,        # ticks along the left edge are off
        labelbottom=False, # labels along the bottom edge are off
        labelleft=False    # labels along the left edge are off
    )
    plt.grid()
    plt.axis([xmin,xmax,ymin,ymax])
    plt.title('Simulation Time = '+str(int(t)))
    plt.savefig(output_filename,format="png")
    plt.close(fig)

    ### evolve the particle's position ###
    for i in range(0,npart):     
        ### take a tentative step ###
        vx[i] += ax[i]*dt
        vy[i] += ay[i]*dt
        pxtemp = px[i] + vx[i]*dt
        pytemp = py[i] + vy[i]*dt

        ### check if collision with x = xmax wall  ###
        ### modify px, py, vx, vy if collision     ###
        if pxtemp+r[i]>xmax:
            vx[i] = -1.0*vx[i]
            pxtemp = xmax-r[i]  ### the collision takes 1 time interval

        ### check if collision with x = xmin wall  ###
        ### modify px, py, vx, vy if collision     ###
        if pxtemp-r[i]<xmin:
            vx[i] = -1.0*vx[i]
            pxtemp = xmin+r[i]  ### the collision takes 1 time interval
        
        ### check if collision with y = ymax wall  ###
        ### modify px, py, vx, vy if collision     ###
        if pytemp+r[i]>ymax:
            vy[i] = -1.0*vy[i]
            pytemp = ymax-r[i]  ### the collision takes 1 time interval
        
        ### check if collision with y = ymin wall  ###
        ### modify px, py, vx, vy if collision     ###
        if pytemp-r[i]<ymin:
            vy[i] = -1.0*vy[i]
            pytemp = ymin+r[i]   ### the collision takes 1 time interval
        
        px[i] = pxtemp
        py[i] = pytemp

        for j in range(0,npart):
            if i != j:
                vx[i],vy[i],vx[j],vy[j],px[i],py[i],px[j],py[j] = \
                    collision(px[i],px[j],py[i],py[j],r[i],r[j],vx[i],vx[j],vy[i],vy[j],m[i],m[j])

    iframe = iframe+1

### stitch together the frames to make a movie ### 
os.system("rm hw05a.mp4")
os.system("ffmpeg -i %d.png -vf scale=432x432 hw05a.mp4")
os.system("rm *.png") 


'''
############## Choice 2 ##############
#!/usr/plocal/bin/python3

### PHYS 5730/7730 spring 2022 ###
### Homework #5(a)             ### 
### JB - instructor's solution ###

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import time

matplotlib.rcParams.update({'font.size': 16})

### animation parameters ###

xmin =  0.0
xmax = 10.0
ymin =  0.0
ymax = 10.0
tmin =  0.0
tmax = 60.0
vmin =  0.0
vmax =  1.0
dt   =  0.05

gravity = False
npart   = 20

### define particle(s) ###

### (nb - ideally particles should move less than 1 radius in dt) ### 

npart = 20

px = np.random.uniform(xmin+1.0,xmax-1.0,npart)
py = np.random.uniform(ymin+1.0,ymax-1.0,npart)

vx = np.random.uniform(-vmax,vmax,npart)
vy = np.random.uniform(-vmax,vmax,npart)

#r  = np.random.uniform(0.1,0.5,npart) # particle radius in sim units
r  = 0.2*np.ones(npart)               # particle radius in sim units  
rp = 65*r                             # particle radius in points

m  = 1.0*np.ones(npart)*r*r*r

colors  = ['b']*npart
colors[0] = 'r'

### turn on gravity? ###

if(gravity):
    g  = -0.5
    ax = np.zeros(npart)
    ay = g*np.ones(npart)
else:
    g  = -0.5
    ax = np.zeros(npart)
    ay = np.zeros(npart)

#######################    
def collide(i,j,dtest):
#######################

#   2D elastic collision

    if(True): 

        M = m[i]+m[j]

        dpxij = px[i] - px[j]
        dpyij = py[i] - py[j]
        dpxji = px[j] - px[i]
        dpyji = py[j] - py[i]

        dvxij = vx[i] - vx[j]
        dvyij = vy[i] - vy[j]
        dvxji = vx[j] - vx[i]
        dvyji = vy[j] - vy[i]

        dvpij = (vx[i]-vx[j])*(px[i]-px[j]) + (vy[i]-vy[j])*(py[i]-py[j]) 
        dvpji = (vx[j]-vx[i])*(px[j]-px[i]) + (vy[j]-vy[i])*(py[j]-py[i]) 
    
        dsq  = dpxij*dpxij + dpyij*dpyij
    
        vxinew = vx[i] - (2.0*m[j]/M)*dvpij*dpxij/dsq
        vyinew = vy[i] - (2.0*m[j]/M)*dvpij*dpyij/dsq

        vxjnew = vx[j] - (2.0*m[i]/M)*dvpji*dpxji/dsq
        vyjnew = vy[j] - (2.0*m[i]/M)*dvpji*dpyji/dsq

        vx[i] = vxinew
        vy[i] = vyinew

        vx[j] = vxjnew
        vy[j] = vyjnew

        ### also tweak positions so particles don't stick ###
        
        ### find unit vector along line joining centers of i, j ###

        vijx = vx[j] - vx[i]
        vijy = vy[j] - vy[i]
        vunx = vijx/np.sqrt(vijx*vijx + vijy*vijy)
        vuny = vijy/np.sqrt(vijx*vijx + vijy*vijy)

        ### relocate particles along this unit vector ###

        px[i] = px[i] - vunx*abs(dtest)
        py[i] = py[i] - vuny*abs(dtest)
        px[j] = px[j] + vunx*abs(dtest)
        py[j] = py[j] + vuny*abs(dtest)
        
    return

### create the current frame, evolve particle positions ### 
    
iframe = 0

for t in np.arange(tmin,tmax,dt):

    ### output the current frame ###
    
    output_filename = '%d.png' % iframe
    
    fig = plt.figure(figsize=(6,6))

    plt.scatter(px,py,c=colors,s=rp*rp)
    
    plt.tick_params(
        axis='both',       # changes apply to both axes
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        right=False,       # ticks along the right edge are off
        left=False,        # ticks along the left edge are off
        labelbottom=False, # labels along the bottom edge are off
        labelleft=False    # labels along the left edge are off
    ) 

    plt.grid()
    plt.axis([xmin,xmax,ymin,ymax])
    plt.title('Simulation Time = '+str(int(t)))
    plt.savefig(output_filename,format="png")
    plt.close(fig)

    ### evolve the particles' positions ###

    for i in range(0,npart):

        ### take a tentative step ###
        
        pxtemp = px[i] + vx[i]*dt + 0.5*ax[i]*dt*dt
        pytemp = py[i] + vy[i]*dt + 0.5*ay[i]*dt*dt

        vxtemp = vx[i] + ax[i]*dt
        vytemp = vy[i] + ay[i]*dt

        ### Check if a collision with a wall occurred. If so, reverse ###
        ### sign of perpendicular component of velocity. Also, insure ###
        ### 'bounce' large enough to avoid wall-stickery.             ###
        
        if(pxtemp+r[i] > xmax):
            vxtemp = -1.0*vxtemp
            pxtemp = 2.0*(xmax-r[i]) - pxtemp
            
        if(pxtemp-r[i] < xmin):
            vxtemp = -1.0*vxtemp
            pxtemp = 2.0*(xmin+r[i]) - pxtemp

        if(pytemp+r[i] > ymax):
            vytemp = -1.0*vytemp
            pytemp = 2.0*(ymax-r[i]) - pytemp
    
        if(pytemp-r[i] < ymin):
            vytemp = -1.0*vytemp 
            pytemp = 2.0*(ymin+r[i]) - pytemp

        px[i] = pxtemp
        py[i] = pytemp
        vx[i] = vxtemp
        vy[i] = vytemp

        ### Check if collision with another particle occurred. ###
        ### If so, call collision function above.              ###

        for j in range(0,npart):

            if(i==j):
                continue

            distij = np.sqrt((px[i]-px[j])*(px[i]-px[j])
                             + (py[i]-py[j])*(py[i]-py[j]))

            dtest = distij - (r[i]+r[j])
            
            if(dtest < 0):     
                collide(i,j,dtest)
                break

    iframe = iframe+1
    
### stitch together the frames to make a movie ### 
        
os.system("rm hw05a.mp4")
os.system("ffmpeg -i %d.png -vf scale=432x432 hw05a.mp4")
os.system("rm *.png") 

#
'''