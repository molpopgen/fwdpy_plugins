import pyximport
pyximport.install()
import snowdrift as sd
import fwdpy as fp
import fwdpy.fitness as fpw
import numpy as np
import sys
#sys.exit(1)
nregions=[]
sregions=[fp.GaussianS(0,1,1,0.25)]
recregions=[fp.Region(0,1,1)]

N=1000
nlist=np.array([1000]*(10*N),dtype=np.uint32)
n=0.0
s=0.001
r=0.01
w=sd.SpopSnowdrift(-0.1,0.25,0.3,-1.0)
sampler=fp.NothingSampler(1)
pops=fp.SpopVec(1,N)
rng=fp.GSLrng(101)
for i in pops:
    sd.label_diploids(i)
fp.evolve_regions_sampler_fitness(rng,pops,sampler,w,
                                  nlist,
                                  n,s,r,
                                  nregions,sregions,recregions,0)
#print( "here")
for i in pops:
    view=fp.view_diploids(i,[0,1,2,4,5,6])
    for i in view:
        print("f",i['w'])#,i['g'])
