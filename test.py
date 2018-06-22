
import FHNC_3D as FHNC
import numpy as np

FHNC.set_default()
FHNC.nu=2
FHNC.update()
gr=FHNC.solve_EL(10.)

tmp=np.array([FHNC.x,gr])
np.savetxt("result.txt",tmp.transpose(),fmt='%1.5e',header='#r  g(r)')
