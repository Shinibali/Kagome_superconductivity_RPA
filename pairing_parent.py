#!/usr/bin/env python

import numpy as np
import subprocess as sp
import time as t
import sys
import struct
import os

''' change fdid info in unique_bare_interactions line 7 '''
fdid = '/data/shinibali/kagome/oneorb/0p08/'
directory = os.path.dirname(fdid)
if not os.path.exists(directory):
    os.makedirs(directory)

ncont = 6*25 # should be multiple of 6 for kf along a hexagonal ring
## U =   [  -1.0, 0.0 ]
## Vnn = [  0.0, -1.0 ]
## J =   [  0.1, 0.1 ]

U =   [0.0, 0.4, 0.4, 0.8, 0.8, 0.8, 1.2, 1.2, 1.2, 1.2, \
       1.6, 1.6, 1.6, 1.6, 2.0, 2.0, 2.0, 2.0 ]
Vnn = [0.0, 0.0, 0.2, 0.0, 0.2, 0.4, 0.0, 0.2, 0.4, 0.6, \
       0.0, 0.2, 0.4, 0.6, 0.0, 0.2, 0.4, 0.6 ]
J =   [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, \
       0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ]

t0 = t.time()
print('Starting to interpolate bare susceptibility (MATLAB codes).\n')
 
for n in range(0,len(U)):
    print('Starting to evaluate interaction matrices for U_{0}(MATLAB codes).\n').format(n+1)
    uu = U[n]
    jj = J[n]
    vv = Vnn[n]
    inputs = [ncont,uu,jj,vv]
    fun = '-r "pairing_bare_interactions(%s); quit" >/dev/null ' % (str(inputs))
    lmb = ['/usr/local/bin/matlab',' -nodesktop',' -nosplash',' -nodisplay',fun]
    sp.Popen(lmb).wait()

    t4 = t.time()-t0
    print("Finished evaluating interaction matrices for U_{0}. Time elasped: {1:.2f} s, or {2:.2f} hr.").format(n+1,t4,t4/3600.)
    print('Starting calculation of pair vertex in band space using f90 code.\n')

    lmb = ['/usr/bin/gfortran','-O3','-march=native','-ffast-math']
    tus = ['./functions/module.f90','./functions/pv_eval.f90','-llapack','-lblas']
    sp.Popen(lmb+tus).wait()
    sp.Popen(['./a.out',str(fdid),str(n+1)]).wait()

    t4 = t.time()-t0
    print("Finished evaluating pv in band space for U_{0}. Time elasped: {1:.2f} s, or {2:.2f} hr.\n").format(n+1,t4,t4/3600.)

print("Finished all U sets. Collecting and organizing data sets.\n")

fun = '-r "linearized_gap(%s,%s,%s); quit" >/dev/null ' % (str(U),str(J),str(Vnn))
lmb = ['/usr/local/bin/matlab',' -nodesktop',' -nosplash',' -nodisplay',fun]
sp.Popen(lmb).wait()

t4 = t.time()-t0
print("Finished all calculations. Exiting python. Time elasped: {0:.2f} s, or {1:.2f} hr.\n").format(t4,t4/3600.)
