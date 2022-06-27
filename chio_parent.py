#!/usr/bin/env python

import numpy as np
import subprocess as sp
import time as t
import sys
import struct
import os

qxsus = 80          # divisible by 3, q-grid for susceptibility
kxsus = 2*qxsus     # k-grid for sum to evaluate susceptibility
orbitals, Ef = 1,0.08  # chemical potential
MatsuT = 116.05     # this goes as k*MatsuT for Temperature in Kelvin
nproc = 45

t0 = t.time()
fdic = '/data/shinibali/kagome/oneorb/0p08/'
fres = './bareresults/'
fin = './input_files/'
directory = os.path.dirname(fres)
if not os.path.exists(directory):
    os.makedirs(directory)
directory = os.path.dirname(fin)
if not os.path.exists(directory):
    os.makedirs(directory)
directory = os.path.dirname(fdic)
if not os.path.exists(directory):
    os.makedirs(directory)

print('\nStarting matlab engine.\n')
print('Writing input files for: {0}-{1} k-points.\n').format(kxsus,kxsus)

inputs = [kxsus,qxsus,orbitals,Ef,MatsuT] 
fun = '-r "chio_input(%s); quit" >/dev/null ' %( str(inputs) )
lmb = ['/usr/local/bin/matlab',' -nodesktop',' -nosplash',' -nodisplay',fun]
sp.Popen(lmb).wait()

t2 = t.time()-t0
print('Finished writing input files. Time elapsed: {0:.2f} s, or {1:.2f} hr.\n').format(t2,t2/3600.)

totq = int(qxsus*qxsus) # total no. of q-points for bare susceptibility
iq = range(1,totq+1)
gp = totq/nproc
if totq%nproc != 0:
   gp = gp+1
list = [iq[j:j+gp] for j in range(0,len(iq),gp)]

print('Starting child processes for calculating bare susceptibility at each q-point.\n')
print('Distributing workload to {} processors with {} q-points in each.\n').format(len(list),gp)
est_time = 0.000005019 * gp * (kxsus)**2 * (orbitals *3 )**6
print('Estimated Time consumption: {0:.2f} s, or {1:.2f} hr.\n').format(est_time,est_time/3600.)

t6 = t.time()-t0
lmb = ['/usr/bin/gfortran','-O3','-march=native','-ffast-math']
tus = ['./functions/module.f90','./functions/susceptibility.f90','-llapack','-lblas']
sp.Popen(lmb+tus).wait()

p=[]
for n in range(0,len(list)): 
    ik = list[n]
    lmb = ['./a.out',str(ik[0]),str(ik[-1])]
    p.append(sp.Popen(lmb))
for q in p:
    q.wait()

t5 = t.time()-t0
print('All of total {0} batches of q-points finished execution. Time elapsed: {1:.2f} s, or {2:.2f} hr.').format(len(list),t5,t5/3600.)
print('Individual timing for SE parallelization: {0:.2f} s, or {1:.2f} hr.\n').format(t5-t6,(t5-t6)/3600.)
print('Starting to collect and organize bare susceptibility data using f90 code.\n')

lmb = ['/usr/bin/gfortran','-O3','-march=native','-ffast-math','./functions/collect_baresusdata.f90']
sp.Popen(lmb).wait()
sp.Popen('./a.out').wait()

fun = '-r "chio_arrange(%s); quit" >/dev/null ' % (str(qxsus))
lmb = ['/usr/local/bin/matlab',' -nodesktop',' -nosplash',' -nodisplay',fun]
sp.Popen(lmb).wait()

t3 = t.time()-t0
print('All data collected and organized. Time elapsed: {0:.2f} s, or {1:.2f} hr.').format(t3,t3/3600.)
