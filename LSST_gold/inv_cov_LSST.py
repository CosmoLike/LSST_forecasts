#!/usr/bin/python
import sys
sys.path.insert(0, '/home/teifler/Python-2.7.8/')
import math, numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy import linalg as LA
import numpy as np

infile = 'cov/Multi_Probe_LSST_2.600000e+01_1.800000e+04_Rmin10_lmin30_lens-zmin0.15LSST_cov_Ncl20_Ntomo10'
#infile = 'Multi_Probe_LSST_2.600000e+01_1.800000e+04_Rmin10LSST_cov_Ncl20_Ntomo10_Nlens10_noredshifterr_old_bias_scale_g_inv'
nggl = 42 # number of lens-source bins. Since a bin is not used when lens has siginificant overlap with source, it is given in code output
ngcl = 0 # 
nlens = 10 # number of lens bin
nlenscl=0 # 3 if set_cluster_DES, 4 if set_cluster_LSST 
nshear = 55 # shear x shear 
ncl=20 # number of cl shared among galaxy-galaxy lensing, galaxy-galaxy clusetering, and cosmic shear
nclgcl=0
nrich=0

ndata = (nshear+nggl+nlens)*ncl+nlenscl*nrich+nrich*ngcl*nclgcl #??? for WFIRST
n2pt = (nshear+nggl+nlens)*ncl #(55+?+4)*20 for WFIRST
ncluster = nlenscl*nrich #4*7 fir WFIRST
n2ptcl=n2pt+ncluster
nclusterN_WL=ncluster+nrich*ngcl*nclgcl
# ncluster =0
# n2ptcl=0
# nclusterN_WL=0
covfile = np.genfromtxt(infile)
cov = np.zeros((ndata,ndata))

print ndata,n2pt,int(np.max(covfile[:,0])+1)

for i in range(0,covfile.shape[0]):
  cov[int(covfile[i,0]),int(covfile[i,1])] = covfile[i,8]+covfile[i,9]
  cov[int(covfile[i,1]),int(covfile[i,0])] = covfile[i,8]+covfile[i,9]

filename = infile + "_sorted"
np.savetxt(filename, cov)
print "sorted covariance written in:", filename 

# for i in range(n2ptcl,ndata):
#   for j in range(0,nshear*ncl):  
#     cov[i,j]=0.0
#     cov[j,i]=0.0

# for i in range(0,covfile.shape[0]):
#   cov[int(covfile[i,0]),int(covfile[i,1])] += covfile[i,9]
#   cov[int(covfile[i,1]),int(covfile[i,0])] += covfile[i,9]

cor = np.zeros((ndata,ndata))
for i in range(0,ndata):
    for j in range(0,ndata):
       if (cov[i,i]*cov[j,j] >0):
         cor[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])

a = np.sort(LA.eigvals(cor))
print "min+max eigenvalues full cor:"
print np.min(a), np.max(a)
print "neg eigenvalues full cor:"
for i in range(0,ndata):
  if (a[i]< 0.0): print a[i]

b = np.sort(LA.eigvals(cov))
print "min+max eigenvalues full cov:"
print np.min(b), np.max(b)
print "neg eigenvalues full cov:"
for i in range(0,ndata):
  if (b[i]< 0.0): print b[i]


# outfile = infile+"_shear"
# f = open(outfile, "w")
# for i in range(0,nshear*ncl):
#   for j in range(0,nshear*ncl):
#     f.write("%d %d %e\n" %(i,j,cov[i,j]))
# f.close()


# ############### invert shear covariance #################
inv = LA.inv(cov[0:nshear*ncl,0:nshear*ncl])
a = np.sort(LA.eigvals(cov[0:nshear*ncl,0:nshear*ncl]))
print "min+max eigenvalues shear cov:"
print np.min(a), np.max(a)
outfile = infile+"_shear_inv"
f = open(outfile, "w")
for i in range(0,nshear*ncl):
  for j in range(0,nshear*ncl):
    f.write("%d %d %e\n" %(i,j, inv[i,j]))
f.close()

# ############### invert ggl covariance #################
inv = LA.inv(cov[nshear*ncl:(nshear+nggl)*ncl,nshear*ncl:(nshear+nggl)*ncl])
a = np.sort(LA.eigvals(cov[nshear*ncl:(nshear+nggl)*ncl,nshear*ncl:(nshear+nggl)*ncl]))
print "min+max eigenvalues ggl cov:"
print np.min(a), np.max(a)
outfile = infile+"_ggl_inv"
f = open(outfile, "w")
for i in range(0,nggl*ncl):
  for j in range(0,nggl*ncl):
    f.write("%d %d %e\n" %(i,j, inv[i,j]))
f.close()

# ############### invert clustering covariance #################
inv = LA.inv(cov[(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl,(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl])
a = np.sort(LA.eigvals(cov[(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl,(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl]))
print "min+max eigenvalues clustering cov:"
print np.min(a), np.max(a)
outfile = infile+"_clustering_inv"
f = open(outfile, "w")
for i in range(0,nlens*ncl):
  for j in range(0,nlens*ncl):
    f.write("%d %d %e\n" %(i,j, inv[i,j]))
f.close()

# ############### invert ggl + clustering covariance #################
inv = LA.inv(cov[nshear*ncl:(nshear+nggl+nlens)*ncl,nshear*ncl:(nshear+nggl+nlens)*ncl])
a = np.sort(LA.eigvals(cov[nshear*ncl:(nshear+nggl+nlens)*ncl,nshear*ncl:(nshear+nggl+nlens)*ncl]))
print "min+max eigenvalues ggl + clustering cov:"
print np.min(a), np.max(a)
outfile = infile+"_ggl_clustering_inv"
f = open(outfile, "w")
for i in range(0,(nlens+nggl)*ncl):
  for j in range(0,(nlens+nggl)*ncl):
    f.write("%d %d %e\n" %(i,j, inv[i,j]))
f.close()

# ############### invert 2pt covariance #################
a = np.sort(LA.eigvals(cov[0:n2pt,0:n2pt]))
print "min+max eigenvalues 2pt cov:"
print np.min(a), np.max(a)
print "number of negative eigenvalues: %s/%s" % (np.sum(a < 0.), len(a))
inv = LA.inv(cov[0:n2pt,0:n2pt])
outfile = infile+"_2pt_inv" 
f = open(outfile, "w")
for i in range(0,n2pt):
  for j in range(0,n2pt):
    f.write("%d %d %e\n" %( i,j, inv[i,j]))
f.close()


# # ############### invert clusterN covariance #################
# inv = LA.inv(cov[n2pt:n2pt+ncluster,n2pt:n2pt+ncluster])
# a = np.sort(LA.eigvals(cov[n2pt:n2pt+ncluster,n2pt:n2pt+ncluster]))
# print "min+max eigenvalues of clusterN matrix:"
# print np.min(a), np.max(a)
# if (np.min(a)<0):
#   print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile)

# outfile = infile+"_clusterN_inv"
# f = open(outfile, "w")
# for i in range(0,ncluster):
#   for j in range(0,ncluster):
#     f.write("%d %d %e\n" %( i,j, inv[i,j]))
# f.close()


# # ############### invert flull2pt+clusterN+clusterWL covariance #################
# precond = 1.e-7
# for i in range(0,ncluster):
#   cov[n2pt+i,:]*= precond
#   cov[:,n2pt+i]*= precond
# inv = LA.inv(cov)
# a = np.sort(LA.eigvals(cov))
# print "min+max eigenvalues of full 2ptclusterN+clusterWL pre-conditioned matrix:"
# print np.min(a), np.max(a)
# if (np.min(a)<0):
#   print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile)
# for i in range(0,ncluster):
#   inv[n2pt+i,:]*= precond
#   inv[:,n2pt+i]*= precond

# outfile = infile+"_2pt_clusterN_clusterWL_inv"
# f = open(outfile, "w")
# for i in range(0,ndata):
#   for j in range(0,ndata):
#     f.write("%d %d %e\n" %( i,j, inv[i,j]))
# f.close()



# # ############### invert clusterN+clusterWL covariance #################
# inv = LA.inv(cov[n2pt:n2pt+nclusterN_WL,n2pt:n2pt+nclusterN_WL])
# a = np.sort(LA.eigvals(cov[n2pt:n2pt+nclusterN_WL,n2pt:n2pt+nclusterN_WL]))
# print "min+max eigenvalues of clusterN_WL pre-conditioned matrix:"
# print np.min(a), np.max(a)
# if (np.min(a)<0):
#   print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile)

# for i in range(0,ncluster):
#   inv[i,:]*= precond
#   inv[:,i]*= precond

# outfile = infile+"_clusterN_clusterWL_inv"
# f = open(outfile, "w")
# for i in range(0,nclusterN_WL):
#   for j in range(0,nclusterN_WL):
#     f.write("%d %d %e\n" %( i,j, inv[i,j]))
# f.close()

# ############### invert 2pt+clusterN covariance #################
# inv = LA.inv(cov[0:n2pt+ncluster,0:n2pt+ncluster])
# a = np.sort(LA.eigvals(cov[0:n2pt+ncluster,0:n2pt+ncluster]))
# print "min+max eigenvalues of 2pt+N matrix:"
# print np.min(a), np.max(a)
# if (np.min(a)<0):
#   print "WARNING  WARNING: %s is not positive definite! WARNING!" % (infile)

# for i in range(0,ncluster):
#   inv[n2pt+i,:]*= precond
#   inv[:,n2pt+i]*= precond
# outfile = infile+"_2pt_clusterN_inv"
# f = open(outfile, "w")
# for i in range(0,n2pt+ncluster):
#   for j in range(0,n2pt+ncluster):
#     f.write("%d %d %e\n" %( i,j, inv[i,j]))
# f.close()


# cor = np.zeros((ndata,ndata))
# for i in range(0,ndata):
#     for j in range(0,ndata):
#       if (cov[i,i]*cov[j,j] >0):
#         cor[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])

plt.figure()
#plt.imshow(np.log10(cov[0:1500,2000:]), interpolation="nearest",vmin=-25, vmax=-10)
plt.imshow(np.log10(cov[1200:1300,0:100]), interpolation="nearest",vmin=-25, vmax=-10)
#plt.imshow(cor[n2ptcl:n2ptcl+200,300:nshear*ncl], interpolation="nearest",vmax=0.5)
plt.colorbar()
plt.show()
