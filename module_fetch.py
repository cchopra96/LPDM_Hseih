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

def calcfetch(zm,z0,L):
    """
    Calculate the fetch x dependent on instrument height, roughness length, and Monin-Obukhov length
    
    Parameters
    ----------
    
    zm: single element (float)
      instrument height (m)
      
    z0: single element (float) 
      roughness length (m)
      
    L: single element (float)
      Monin-Obukhov length (m)
      
    Returns
    -------
    
    x: single element (float) 
      fetch (m)
      
    other constants like nu, nd, N, k, ustar, x, z
    """
    import numpy as np
    
    ustar=0.12
    N=1000
    k=0.4#von-Karman constant
    xi=0#initial x coordinate
    zi=z0#initial z coordinate
    nu=0
    nd=0
    iter1=56000

    
    
    #initial conds
    x=np.zeros((N,iter1))
    z=np.zeros((N,iter1))
    
    U=np.zeros((N,iter1))
    sig_w=np.zeros((N,iter1))
    dss_dz=np.zeros((N,iter1))
    w=np.zeros((N,iter1))
    a=np.zeros((N,iter1))
    b=np.zeros((N,iter1))
    dw=np.zeros((N,iter1))
    
    x[:,0]=xi
    z[:,0]=zi
    if L >= 0:
        zhi=-5*zm/L
        sig_wi=1.25*ustar
        phi_hi=1+5*zm/L
    
    
    
    elif L < 0:
        y=(1-16*zm/L)**0.25
        zhi=2*np.log((y+1)/2)+np.log((y**2+1)/2)+2*np.arctan((1-y)/(1+y))
        sig_wi=1.25*ustar*(1-3*zm/L)**(1/3)
        phi_hi=0.32*(0.037-zm/L)**(-1/3)
        
        
    tL=(k*zm*ustar)/(phi_hi*sig_wi**2)
    dt=0.1*tL
    ss=1
    i=0   
    
    
    
    while (nu-nd) < 0.9*N:
  
        
        for j in range(N):
       
            U[j,i+1]=np.abs((ustar/k)*(np.log(z[j,i]/z0)-zhi))
            x[j,i+1]=x[j,i]+U[j,i+1]*dt
            
            if L>=0:
                sig_w[j,i]=1.25*ustar
                dss_dz[j,i]=0#second derivative of sig_w_sq;
            elif L<0:
                sig_w[j,i]=1.25*ustar*(1-3*z[j,i]/L)**(1/3)
                dss_dz[j,i]=-3.1250*((ustar**2)/L)*(1-3*z[j,i]/L)**(-1/3)
            
            a[j,i]=-(w[j,i]/tL)+(0.5*dss_dz[j,i]*(1+(w[j,i]/sig_w[j,i])**2))
            b[j,i]=((2*sig_w[j,i]**2)/tL)**0.5
            
            dw[j,i]=a[j,i]*dt+b[j,i]*np.random.normal(0,np.sqrt(dt))
            w[j,i+1]=w[j,i]+dw[j,i]
            z[j,i+1]=z[j,i]+ss*w[j,i]*dt
            
            if (z[j,i]<zm) and (z[j,i+1]>=zm):
                nu=nu+1
            elif (z[j,i]>=zm) and (z[j,i+1]<zm):
                nd=nd+1
                
            if (z[j,i]>=z0) and (z[j,i+1]<z0) and (w[j,i]<0):
                w[j,i]=-w[j,i]
                z[j,i+1]=z[j,i]+ss*w[j,i]*dt
            
            
        i=i+1
        
    
    x[x==0]=np.nan
    z[z==0]=np.nan
    x_fetch=x[j,i-1]
    return x_fetch,nu,nd,N,k,ustar,x,z
