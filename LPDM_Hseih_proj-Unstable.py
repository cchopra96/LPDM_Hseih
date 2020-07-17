# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import module_fetch#library with function for generating fetch
from module_fetch import calcfetch
import numpy as np
import matplotlib.pyplot as plt


# +
#pick case for stability (1 for unstable, 2 for neutral, 3 for stable)
caseA=1 #unstable


zm_s1=np.linspace(2,20,10)
z0_s1=np.linspace(0.01,0.1,10)
#L_s1=([-0.1, -0.08, -0.05, -0.02, -0.01, 1, 5, 10, 20, 50])
L_s1=([-0.1, -0.08, -0.05, -0.01, 30, 70, 100, 150, 200, 250])

l_zm=len(zm_s1)
l_z0=len(z0_s1)
l_L=len(L_s1)
# -

arr_data=[(i,j,k) for i in zm_s1 for j in z0_s1 for k in L_s1]
arr_data1=np.array(arr_data)

zu_over_L_s1=(arr_data1[:,0]/arr_data1[:,2])*(np.log(arr_data1[:,0]/arr_data1[:,1])-1+arr_data1[:,1]/arr_data1[:,0])

np.shape(zu_over_L_s1)

data_unstable=arr_data1[(zu_over_L_s1<-0.04),:]
data_stable=arr_data1[(zu_over_L_s1>0.04),:]
data_neutral=arr_data1[(zu_over_L_s1>-0.04)&(zu_over_L_s1<0.04),:]

print(np.shape(data_unstable))
print(np.shape(data_stable))
print(np.shape(data_neutral))

if caseA==1:
    data=data_unstable
elif caseA==2:
    data=data_neutral
elif caseA==3:
    data=data_stable

# +
kkk,kkk1=data.shape

x_end=np.zeros(kkk)
kkk1=kkk-1

# +

for g in range(kkk):

    z0=data[g,1]
    zm=data[g,0]
    L=data[g,2]

    x_end[g],nu,nd,N,k,ustar,x,z=calcfetch(zm,z0,L)
   
    #the following lines ay be uncommented to plot the ensemble particle trajectories
    #for iii in range(N):
     #   plt.plot(x[iii,:],z[iii,:])
     #   plt.xlabel('x') 
     #   plt.ylabel('z') 
        
    #plt.show()
   
    #(nu-nd)/N
# -


X=np.log((data[:,0]/np.abs(data[:,2]))*(np.log(data[:,0]/data[:,1])-1+data[:,1]/data[:,0]))
Y=np.log(x_end/np.abs(data[:,2]))

plt.scatter(X,Y)
plt.xlabel('log($z_u$/|L|)') 
plt.ylabel('log(x/|L|)') 

p=np.polyfit(X,Y,1);
p

# +
par_P=p[0];
par_D=-(k**2)*np.log(0.9)*np.exp(p[1]);

print([par_P,par_D])
# -

#run this section only for unstable case
if caseA==1:
    dat=[X,Y]
    d1=len(X)
    
    
    for dd in range(d1):
        if X[dd]>6.9 and Y[dd]<6.4:
            X[dd]=np.nan;
            Y[dd]=np.nan;
    X1=X
    Y1=Y 
    
    
    idx = np.isfinite(X1) & np.isfinite(Y1)
    p2 = np.polyfit(X1[idx], Y1[idx], 1)
    
    plt.scatter(X1,Y1)
    plt.xlabel('log($z_u$/|L|)')
    plt.ylabel('log(x/|L|)') 
    
    par_P1=p2[0];
    par_D1=-(k**2)*np.log(0.9)*np.exp(p2[1]);
    
    print([par_D1,par_P1])




