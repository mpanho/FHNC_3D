
import FHNC_3D as FHNC
import numpy as np

rs=input('Enter rs vaule:')
FHNC.set_default()

pol=input('for spin polarized calculation enter [y]:')
if (pol=='y'):
  FHNC.nu=1
  FHNC.update()

gr=FHNC.solve_EL(rs)

tmp=np.array([FHNC.x,gr])
np.savetxt("result.txt",tmp.transpose(),fmt='%1.5e',header='#r/rs/a0  g(r) ')
print('Output written to: result.txt')
