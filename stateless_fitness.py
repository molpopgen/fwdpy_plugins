from __future__ import print_function
import pyximport
pyximport.install()
import stateless_fitness_models as sfm
import numpy as np
import fwdpy as fp


n = [fp.Region(0,1,1)]
s = [fp.ExpS(0,1,1,-0.25)]
r = n

def evolve_some_pops(w):
    print ("starting!")
    N=1000

    rng = fp.GSLrng(101)
    pops = fp.SpopVec(10,N)
    nlist = np.array([N]*(10*N),dtype=np.uint32)
    sampler = fp.NothingSampler(len(pops))
    fp.evolve_regions_sampler_fitness(rng,pops,sampler,w,nlist,0.,0.0025,0.01,n,s,r,0)

    dips1 = fp.view_diploids(pops,[0,1,2,3,4,5])

    for i in dips1[0]:
        print (i['w'])

for w in [sfm.CustomAdditive1(),sfm.CustomAdditive2(),sfm.CustomAdditive3(),sfm.CustomAdditive4()]:
    evolve_some_pops(w)
    







































