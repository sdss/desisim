
from astropy.io import fits as pyfits 
import numpy as np
import matplotlib.pyplot as plt
hdulist=pyfits.open('sorted_true_elgs.fits')
zarray=hdulist[1].data['ZBTRUEZ']
diff=1.
number_templates=len(zarray)
abs_diff=np.ones(number_templates-1)
diff=1.
for i in range(len(zarray)-1):
    abs_diff[i]=np.abs(zarray[i]-zarray[i+1])
    if abs_diff[i]==0:
        print ("zero diff")
        print i, zarray[i],zarray[i+1],diff
        print hdulist[1].data[i]
        print hdulist[1].data[i+1]
    elif(abs_diff[i]<diff):
        diff=abs_diff[i]
        print i,diff
print diff 
print(" min %f max %f "%(np.amin(zarray),np.amax(zarray)))

#plt.hist(abs_diff,100,(0.,0.001))
#plt.show()   
# make binned list
bin_delta=0.5e-3
number_bins=int((np.amax(zarray)-np.amin(zarray))/bin_delta)

A=[]
x=np.amin(zarray)
B=[]
for i in range(number_templates):
    if zarray[i]<x+bin_delta:
        B.append(i)

    else:
        if((len(B))>0):
            print B
        else:
            B.append(i)
        A.append(B)
        B=[]
        print A
        x=x+bin_delta
for i in range(number_bins):
    print i,A[i]

for i in range(number_bins):
    print i,len(A[i])  
n_per_bin=np.zeros(number_bins) 
for i in range(number_bins):
    n_per_bin[i]=len(A[i])
#plt.hist(n_per_bin,20)
#plt.show()
for k in range(100):
    z=raw_input("Input trial z")

    bin_find=int((z-np.amin(zarray)/bin_delta))
    
    choices=A[bin_find]
    if len(choices)==1:
        z_found=A[0]
    else:
        for j in range(len(choices)):
            if z<choices[j]:
                z_found=choices[j]
            else:
                break
    print z, z_found
