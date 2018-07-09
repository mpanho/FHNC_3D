
import FHNC_3D as FHNC
import numpy as np

tmp=input('Enter rs vaule:')
rs=float(tmp)
FHNC.set_default()

pol=input('for spin polarized calculation enter [y]:')
if (pol=='y'):
  FHNC.nu=1
  FHNC.update()

gr=FHNC.solve_ladder(rs)
Sq=1+FHNC.FFT(gr-1)

kallio=input('for an calculation with the Kallio method enter [y]:')
if (kallio=='y'):
    gr_K=FHNC.solve_Kallio(rs)
    Sq_K=1+FHNC.FFT(gr_K-1)

    tmp=np.array([FHNC.x,gr,gr_K])
    tmp2=np.array([FHNC.k,Sq,Sq_K])
else:
    tmp=np.array([FHNC.x,gr])
    tmp2=np.array([FHNC.k,Sq])


sFHNC=input('for an sFHNC calculation enter [y]:')
if (sFHNC=='y'):
    gr_s=FHNC.solve_sFHNC(rs)
    Sq_s=1+FHNC.FFT(gr_s-1)
    tmp=np.append(tmp,[gr_s],axis=0)
    tmp2=np.append(tmp2,[Sq_s],axis=0)

np.savetxt("result.txt",tmp.transpose(),fmt='%1.5e',header='#r/rs/a0  g(r) ')
print('Output written to: result.txt')
#np.savetxt("Sq.txt",tmp2.transpose(),fmt='%1.5e',header='#q  S(q)')
