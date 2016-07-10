from __future__ import print_function
import pyximport
pyximport.install()
import change_selection_coefficient as csp
import custom_temporal_sampler1 as cts1

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
    print (x)

#Get the list of mutations in each population
views = [pd.DataFrame(i) for i in fp.view_mutations(pops)]

for i in views:
    print (i[i.neutral==False])

#Now, lets use our custom temporal sampler
sampler=cts1.SelectedSFSSampler(len(pops))

#Evolve these pops for another 100 generations,
#tracking joint frequency,s of all selected mutations
fp.evolve_regions_sampler(rng,pops,sampler,
                          nlist[:100], 
                          mutrate_neutral,
                          mutrate_sel,
                          recrate,
                          nregions,
                          sregions,
                          recregions,
                          1)

def coerce2DF(x):
    """
    The data coming out of our sampler are a bit untidy.
    Let's tidy them up and put them in a pandas.DataFrame.
    We'll get the count,s,generation for each generation
    That our selected mutation existed
    """
    generations=[i[0] for i in x if len(i[1])]
    tuples=[i[1][0] for i in x if len(i[1])]
    rv=pd.DataFrame(tuples,columns=['count','s'])
    rv['generation']=generations
    return rv

sfs=[coerce2DF(i) for i in sampler.get()]

for sfsi in sfs:
    print (sfsi)
