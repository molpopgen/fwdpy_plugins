import pyximport
pyximport.install()
import change_selection_coefficient as csp

import fwdpy as fp
import numpy as np
import pandas as pd

nregions=[fp.Region(0,1,1)]
sregions=[]
recregions=nregions

N=1000
nlist=np.array([N]*10*N,dtype=np.uint32)

theta=100.0
rho=100.0

mutrate_neutral = theta/float(4*N)
mutrate_sel = 0.0
recrate = rho/float(4*N)

rng = fp.GSLrng(101)

#Evolve 20 populations to mutation-drift equilibrium
pops = fp.evolve_regions(rng,20,
                         nlist[0],
                         nlist,
                         mutrate_neutral,
                         mutrate_sel,
                         recrate,
                         nregions,
                         sregions,
                         recregions)

#Change a mutation's selection coefficient
for i in pops:
    x=csp.change_selection_coeff(i,0.5,10,0.1,0.25)
    #print out count and position of variant actually changed
    print x

#Get the list of mutations in each population
views = [pd.DataFrame(i) for i in fp.view_mutations(pops)]

for i in views:
    print i[i.neutral==False]
    
