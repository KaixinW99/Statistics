### import package ###
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from array import *
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import scipy.sparse as SciSpa
from numba import jit
import time
import os
import sys

### constants ###
x0 = -3.
k = 50.
sig = 0.25
Cons = np.sqrt((1.0/sig/np.sqrt(np.pi)))  # normalization constant

### define the spacetime ###
xmin = -8.0   ### enlarge the matrix ###
xmax =  8.0   ### enlarge the matrix ###
xstp = 0.001

tmin =  0.0
tmax =  0.2
tstp = 0.001

x = np.arange(xmin,xmax,xstp)
t = np.arange(tmin,tmax,tstp)
Nx = len(x)
Nt = len(t)


### potential energy ###
V = np.zeros(Nx)
V[(x>=-0.5) & (x<=0.5)]=1000
#print(V)
V_mat=SciSpa.diags(V,shape=(Nx,Nx))
#print (V_mat)
### starting wave function ###
psi_init = Cons*np.exp(-(x-x0)**2/2/sig**2)*np.exp(1j*k*(x-x0))
### check the normalization ###
print('Integrated probability: ',np.sum(np.absolute(psi_init)**2*xstp))
### the coefficient of each wavefunction ###
coeff_mat=SciSpa.diags([1,-2,1],[-1,0,1],shape=(Nx,Nx))

### periodic boundary ###
#coeff_mat=SciSpa.csr_matrix(coeff)
#coeff_mat[0,-1]=1.
#coeff_mat[-1,0]=1.

### reflection boundary ###
#coeff_mat=SciSpa.csr_matrix(coeff)
#coeff_mat[0,0] = 2.
#coeff_mat[0,1] = -5.
#coeff_mat[0,2] = 4.
#coeff_mat[0,3] = -1.
#coeff_mat[-1,-4] = -1.
#coeff_mat[-1,-3] = 4.
#coeff_mat[-1,-2] = -5.
#coeff_mat[-1,-1] = 2

### pick larger size ###
print('Coefficient Matrix with size %dx%d'%(Nx,Nx))
print(coeff_mat.toarray())

#@jit(parallel=True)
### cannot use jit to speed up the calculation ###
def SE(t,y,xstp,coeff_mat,V_mat):
    return 1j/2/xstp**2*coeff_mat*y-1j*V_mat*y

#sol = odeint(func=SE,y0=psi_init,t=t,args=(xstp,coeff_mat,V),tfirst=True)
### odeint cannot solve the problems with complex number ###
print('The start to solve PDE')
tstart = time.time()
sol = solve_ivp(fun=SE,t_span=[tmin,tmax],y0=psi_init,t_eval=t,args=(xstp,coeff_mat,V_mat),method='RK23')
tstop = time.time()
print("time to solve  = ",tstop-tstart,"seconds")

print('The start to plot')
tstart = time.time()
### making the frames of the animation ###
for iframe in np.arange(0,len(sol.t)):

    output_filename = '%d.png' % iframe
    titlestring = "time = %5.3f" % t[iframe]

    u = sol.y[:,iframe].real
    v = sol.y[:,iframe].imag

    fig,ax=plt.subplots(1,2,figsize=(20,8))
    ax[0].plot(x,u,label='real')
    ax[1].plot(x,v,label='imag')
    ax[0].set_ylabel('$u(t,x)$')
    ax[1].set_ylabel('$v(t,x)$')
    for axes in ax.flat:
        axes.set_xlabel('$x$')
        axes.legend(loc=1)
        axes.set(xlim=[-5,5],ylim=[-1,1])
        axes.set_title(titlestring)
    plt.savefig(output_filename,format="png")
    plt.close()
tstop = time.time()
print("time to plot  = ",tstop-tstart,"seconds")

### stitching frames together to make animation ###
os.system("rm schro1d.mp4")
os.system("ffmpeg -i %d.png -vf scale=800x800 schro1d.mp4")
os.system("rm *.png") 