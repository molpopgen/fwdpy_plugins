#A simple plugin to make sure that fwdpy typedefs, etc., are cimport-able
#from libcpp
from libcpp.vector cimport vector
from libcpp.utility cimport pair

#From fwdpp
from fwdpy.fwdpp cimport popgenmut

#From fwdpy
from fwdpy.fwdpy cimport gamete_t
from fwdpy.fwdpy cimport gcont_t
from fwdpy.fwdpy cimport ucont_t
from fwdpy.fwdpy cimport mcont_t
from fwdpy.fwdpy cimport diploid_t
from fwdpy.fwdpy cimport dipvector_t
from fwdpy.fwdpy cimport singlepop_t
from fwdpy.fwdpy cimport Spop

#A very boring plugin indeed!
cdef void my_plugin_function(const singlepop_t * pop):
    pass

#This will be the function that your plugin exposes
#to Python:
def foo(Spop p):
    my_plugin_function(p.pop.get())

cdef unsigned count_muts(const singlepop_t * pop) nogil:
    cdef size_t i=0
    cdef size_t n = pop.mcounts.size()
    cdef unsigned twoN = 2*pop.popsize() #This is a member function that returns pop.N
    cdef unsigned extant=0
    for i in range(n):
        #Check that mutation is not extinct and not fixed	
        if pop.mcounts[i] > 0 and pop.mcounts[i] < twoN:
            extant+=1
    #return out count
    return extant
   
cdef unsigned count_neutral_muts(const singlepop_t * pop) nogil:
   cdef size_t i=0
   cdef size_t n = pop.mcounts.size()
   cdef unsigned twoN = 2*pop.popsize() #This is a member function that returns pop.N
   cdef unsigned extant=0
   for i in range(n):
       #Check that mutation is not extinct and not fixed	
       if pop.mcounts[i] > 0 and pop.mcounts[i] < twoN and pop.mutations[i].neutral is True:
            extant+=1
   #return out count
   return extant

#This is a helper function.  It will count the number of segregating mutations
#in each gamete
cdef int count_mutations(const vector[size_t] & keys,const ucont_t & mcounts,const unsigned twoN) nogil:
    cdef size_t i=0
    cdef size_t n=keys.size()
    cdef int rv = 0
    for i in range(n):
        if mcounts[keys[i]] < twoN:
            rv+=1
    return rv
		
cdef vector[pair[int,int]] mutations_per_gamete(const singlepop_t * pop) nogil:
    cdef vector[pair[int,int]] rv
    cdef size_t i = 0
    cdef size_t n = pop.gametes.size()
    cdef unsigned twoN = 2*pop.popsize()
    cdef int neutral,selected
    #Now, we go through every gamete and:
    #1. Check that it is not extinct
    #2. Go over every mutation in each gamete and make sure that it is not fixed.
    #   We do not need to check that each mutation in each gamete has a nonzero count.
    #   fwdpp ensures that an extant gamete contains extant mutations.
    for i in range(n):
        if pop.gametes[i].n > 0: #gamete is not extinct
             neutral = count_mutations(pop.gametes[i].mutations,pop.mcounts,twoN)
             selected = count_mutations(pop.gametes[i].smutations,pop.mcounts,twoN)
             rv.push_back(pair[int,int](neutral,selected))
    return rv

cdef pair[int,int] count_mutations_diploid(const diploid_t & dip, const gcont_t & gametes, const ucont_t & mcounts, const unsigned twoN) nogil:
    cdef int neutral = count_mutations(gametes[dip.first].mutations,mcounts,twoN)
    cdef int selected = count_mutations(gametes[dip.first].smutations,mcounts,twoN)
    return pair[int,int](neutral,selected)

cdef vector[pair[int,int]] mutations_per_diploid(const singlepop_t * pop) nogil:
    cdef vector[pair[int,int]] rv
    cdef size_t i=0
    cdef size_t n=pop.diploids.size()
    cdef unsigned twoN = 2*n
    cdef pair[int,int] temp
    #Now, go through every diploid:
    for i in range(n):
        temp = count_mutations_diploid(pop.diploids[i],pop.gametes,pop.mcounts,twoN)
        rv.push_back(temp)
    return rv

#Population mean fitness under a multiplicative model
import numpy as np;
from cython.view cimport array as cvarray
from fwdpy.fwdpp cimport multiplicative_diploid

cdef void wbar_multiplicative_details(const singlepop_t * pop, double[:] w, const double scaling) nogil:
    cdef multiplicative_diploid wfxn
    cdef size_t i=0, n=pop.diploids.size()
    for i in range(n):
        #Here is the trick.  wfxn is a C++ class, but it is also a function!
        #Further, it is a template function.  Cython is not willing to just let
        #the C++ compiler figure out the types here, so we have to explicitly use typecasts,
        #which is what the <foo>bar is: type cast a bar to a foo.  This has NO RUNTIME PENALTY!!!
        #Yes, we also have to cast the scaling parameter, even though it is not a template parameter.
        w[i] = wfxn(<diploid_t>pop.diploids[i],<gcont_t>pop.gametes,<mcont_t>pop.mutations,<double>scaling)

def wbar_mutiplicative(Spop p, const double scaling):
    w=np.array(p.popsize(),dtype=np.float64)
    wbar_multiplicative_details(p.pop.get(),w[:],scaling)
    return w.mean()
