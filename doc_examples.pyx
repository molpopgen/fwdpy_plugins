#A simple plugin to make sure that fwdpy typedefs, etc., are cimport-able

#From fwdpp
from fwdpy.fwdpp cimport popgenmut

#From fwdpy
from fwdpy.fwdpy cimport gamete_t
from fwdpy.fwdpy cimport gcont_t
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
