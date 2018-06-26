
import numpy as np
from scipy.fftpack import fft, ifft,dst,diff

def set_default():
    global n,xmax,x,dx, k, dk, pi, nu, r0kF, imax
    n=3000
    xmax=50.
    dx=xmax/n
    x=np.arange(0,xmax,xmax/n)+dx/2
    
    pi=np.pi
    dk=pi/xmax#*(n)/(n+1)
    k=np.arange(0,dk*n,dk)  + dk/2

    nu=2  # degeneracy factor, for spin polarized calculation set to 1
    r0kF=(9*pi/2/nu)**(1/3.)
    imax=400
    
def update():
    global n,xmax,x,dx, k, dk, nu, r0kF
    print('new values for n=', n,' xmax=',xmax, ' nu=',nu)
    dx=xmax/n
    x=np.arange(0,xmax,xmax/n)+dx/2
    dk=pi/xmax#*(n)/(n+1)
    k=np.arange(0,dk*n,dk)  + dk/2
    r0kF=(9*pi/2/nu)**(1/3.)


    
def FFT(inp):
    ln=np.arange(n*2)
    inpp=x*inp
    inlong=np.exp(-1j*pi*ln/n/2)*np.concatenate( (inpp,np.flip(inpp,0)) ,0)
    res=np.imag(np.exp(-1j*pi*(1./n/4 + ln/n/2)) * fft(inlong))[0:n]/k
    return -res*dx*3/2#,inlong

def IFFT(inp):
    ln=np.arange(n*2)
    inpp=k*inp
    inlong=np.exp(1j*pi*ln/n/2)*np.concatenate( (inpp,np.flip(inpp,0)) ,0)
    res=np.imag(np.exp(1j*pi*(1./n/4 + ln/n/2)) * ifft(inlong))[0:n]/x
    return n*2*res/xmax/3#*n/(n+1)#,inlong
    
def Diff(inp):
    #res=np.ediff1d(inp,to_end=0)
    res=np.zeros(n)
    res[0]=inp[1]-inp[0]
    for i in range(n-2):
        res[i+1]=inp[i+2]-inp[i]
    
    return res/dx/2
    
def lF(x):
    res=3*(np.sin(x)/x**2-np.cos(x)/x)/x
    return res
    
def SF(x):
    "Free static structure function"
    res= np.heaviside(x-2,0.5)
    res=res +(3*x/4-x**3/16)*np.heaviside(2-x,0.5)
    return res
    
def multip(gin,gold):
    res=np.exp(1*(gin-gold))
    i=0
    while (gin[i]<0.00 and i<n):
        res[i]=np.exp(2*(gin[i]-gold[i]))
        i=i+1
        
    return res
    
def solve_EL(rs):
    "Do the FHNC calculation"
    global g0r,rslist,imax
    print('rs=',rs)
    print('nu=',nu)
    print('If nu=1, spin polarized calculation; if nu=2 paramagnetic calculation')
    print('n=',n)
    print('xmax=',xmax)
    vcoulrs=1*2/x #*np.exp(-x/0.3/xmax)
    gF=1-lF(r0kF*x)**2/nu
    SF0=1+FFT(gF-1)   #SF(k2/r0kF)
    rs_start=.5
    Vphk=FFT(vcoulrs/rs_start/10)
    S0k=SF0 #1/np.sqrt(1 + 2*rs_start**2/k**2 * Vphk)
    g0old=gF
    lap_gFrs= 2/(x**2*np.sqrt(gF)) * Diff(x**2*Diff(np.sqrt(gF)))
    drs=0.5
    tmp=int(rs*2)+2
    rslist=np.linspace(rs_start,rs,tmp)
    j=0
    
    accuracy=4e-4
    hist=[]
    gall=[]
    print('Calculation has started ...')
    for rsi in rslist:
        j=j+1
        i=0
        residual=0.1
        if (rsi==rs):
            accuracy=4e-5
            
        while (residual>accuracy and i<imax):
            i=i+1
            if (i==imax):
                print("not converged at rs=", rsi)
                print("Note: accuracy might be ok for your needs, otherwiese increase imax")
            S0kold=S0k
            wI=-k**2/2/rsi**2 * (1/SF0-1/S0k)**2 * (2*S0k/SF0 +1)
            wIr=IFFT(wI)
            wIB=-k**2/2/rsi**2 * (1-1/S0k)**2 * (2*S0k +1)
            wIBr=IFFT(wIB)
            g00=1+IFFT(S0k-1)
            g0r=g0old*multip(g00,g0old)  #np.exp(1*(g00-g0old))
            g0old=g0r
            Vph= (g0r-0)*vcoulrs/rsi + g0r*(lap_gFrs/rsi**2 +wIr) + (-1)*wIBr + 2/rsi**2*np.abs(Diff(np.sqrt(g0r)))**2
            Vphk=FFT(Vph) + 0*6/k**2/rsi
            if any(np.less_equal(np.real(1 + 2*rsi**2/k**2 * Vphk), np.zeros(n))):
                print("instability at rs=",rsi)
                #break
            S0k=1/np.sqrt(np.abs(1 + 2*rsi**2/k**2 * Vphk))
            g00=1+IFFT(S0k-1)
            residual=np.sum(np.abs(g00-g0r))*dx
            hist.append(residual)
    
        gall.append(g00[0])
    print('...Finished')
    print('residual=',residual, 'iterations=',i)
    return g0r
    
    
