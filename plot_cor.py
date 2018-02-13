#!/usr/bin/python
import sys
import math, numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpig
#from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import linalg as LA
import numpy as np



covfile = './cov_des_y1_i3_unblind_final_v0/COV2_Y1_i3_Ntheta20_Nsource4_Nlens5'

nsource = 4
nshear = nsource*(nsource+1)/2
ntheta = 20
nlens = 5

data = np.genfromtxt(covfile)
ndata = int(np.max(data[:,0]))+1

# additive = np.array([[8.31410482454e-05, 0.000434621003867],
# [0.000164711493213, 0.000455420386097 ],
# [0.000382276276558, -0.000201889775528 ],
# [0.000886358298148, -0.000167289960534 ]])

# additive = np.array([[-0.000490981233482, 0.000371687952072],
# [-8.41019265233e-05, 0.0003731963675],
# [0.000141220372976, -9.89782711765e-05],
# [-0.000123091133757, 0.000621798383881]])

#additive = np.array([[0.000296686456896, 0.000180276352553]])

#additive = np.array([[-0.000151275056988, 0.0002692312903]])
  
# num=0
# xi_sys = np.zeros(nshear)
# for i in range(0,nsource):
#   for j in range(i,nsource):
#     xi_sys[num]=additive[i,0]*additive[j,0]+additive[i,1]*additive[j,1]
#     print i,j,num,xi_sys[num]
#     num += 1
    

# cov_sys= np.zeros((ndata,ndata))
# cov_sys[:,:] = 0.
# for i in range(0,nshear):
#   for j in range(i,nshear):
#     cov_sys[i*20:(i+1)*20,j*20:(j+1)*20]=xi_sys[i]*xi_sys[j]
    
print ndata

ndata_min = int(np.min(data[:,0]))
cov = np.zeros((ndata,ndata))
cov[:,:] = 0.
for i in range(0,data.shape[0]):
  cov[int(data[i,0]),int(data[i,1])] =data[i,8] +data[i,9]
  cov[int(data[i,1]),int(data[i,0])] =data[i,8] +data[i,9]
  
outfile = covfile+".txt"

f = open(outfile, "w")
cor = np.zeros((ndata,ndata))
#cor2 = np.zeros((ndata,ndata))
matrix = np.zeros((ndata,ndata))
for i in range(0,ndata):
  print i, cov[i,i]
  for j in range(0,ndata):
    f.write("%d %d %e\n" %(i,j, cov[i,j]))
    cor[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])
 #   cor2[i,j] = cov2[i,j]/math.sqrt(cov2[i,i]*cov2[j,j])
    # if(i >= j):
    #   matrix[i,j]=cor[i,j]
    # if(i < j):
    #   matrix[i,j]=cor2[i,j]

f.close()

labels = (r'$\xi_+\left(\theta,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',r'$\xi_-\left(\theta,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',r'$\gamma_t\left(\theta;z_\mathrm{l},z_\mathrm{s}\right)$',r'$w\left(\theta;z_\mathrm{l}\right)$')

a = np.sort(LA.eigvals(cor))
b = np.sort(LA.eigvals(cov))
print "min+max eigenvalues cov:"
print b
print np.min(a), np.max(a)
ind = np.where(a < 0)
if (np.min(a) <0): 
	print a[ind]

a = np.sort(LA.eigvals(cor[0:200,0:200]))
print "min+max eigenvalues cov shear++:"
print np.min(a), np.max(a)
ind = np.where(a < 0)
if (np.min(a) <0): 
  print a[ind]

a = np.sort(LA.eigvals(cor[0:400,0:400]))
print "min+max eigenvalues cov shear:"
print np.min(a), np.max(a)
ind = np.where(a < 0)
if (np.min(a) <0): 
	print a[ind]

ticks = np.zeros(5)
tickx = np.zeros(4)
ticks[1] = nshear*ntheta
ticks[2] = 2*nshear*ntheta
ticks[3] = ndata - nlens*ntheta
ticks[4] = ndata
fs= 15
for i in range(0,4):
  tickx[i] = 0.5*(ticks[i]+ticks[i+1])
  plt.plot([ticks[i]-0.5,ticks[i]-0.5],[-.5,ndata-0.5],linestyle ='-',color = 'k')
  plt.plot([-.5,ndata-0.5],[ticks[i]-0.5,ticks[i]-0.5],linestyle ='-',color = 'k')

plt.subplot(1, 1, 1)
ax = plt.gca()
im = ax.imshow(cor[ndata_min:,ndata_min:], interpolation="nearest",origin='lower')
plt.xticks(tickx, labels,fontsize=fs)
plt.yticks(tickx-0.5, labels,fontsize=fs)
plt.tick_params(axis = 'x',length = 0, pad = 15)
plt.tick_params(axis = 'y',length = 0, pad = 5)

plt.colorbar(im)
plt.show()

print ticks
