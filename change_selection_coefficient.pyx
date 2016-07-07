#Lazy: just import everything
from fwdpy.fwdpp cimport *
from fwdpy.fwdpy cimport *
from libc.math cimport *

#Since fwdpy/fwdpp is based on C++,
#Let's use C++ whenever possible
from libcpp.vector cimport vector
from libcpp.limits cimport numeric_limits
from libcpp.utility cimport pair

from cython.operator cimport dereference as deref

#Note: Cython has a weird bug where unsigned
#cannot be used as a type in a C++ pair.
#I still need to report this.
#fwdpy.fwdpy provides uint as a typedef,
#which somehow works.


#This is a Cython function that works on fwdpy/fwdpp
#types.  singlepop_t is defined in include/types.hpp
#in the fwdpy source code.  It publicly inherits from
#fwdpp/sugar/singlepop.hpp, and represents a single-deme,
#contiguout-region population object.
#singlepop_t is exposed to Python via Cython in fwdpy.pxd,
#and is "cimported" into scope above.

#The mechanics of this are straightforward, but do require an
#understanding of what a singlepop_t contains.

#Warning: pop is non-const here, so this is more than enough rope
#to hang oneself with!
cdef pair[uint,double] change_s(singlepop_t * pop,
                                  const double pos,
                                  const unsigned count,
                                  const double new_s,
                                  const double new_h) nogil:
    #Step 1: find the mutation in pop closest to position
    #pos and count "count"
    cdef size_t index = 0
    cdef size_t nmuts = pop.mutations.size()

    #Cython needs concrete objects, so
    #we can't do the more idiomatic
    #auto x = std::numeric_limis<double>::max(), etc.
    cdef numeric_limits[double] dlimits
    cdef numeric_limits[unsigned] ulimits
    cdef double distance = dlimits.max()
    cdef unsigned udistance = ulimits.max()

    #Return now if population is momomorphic
    if nmuts == 0:
        return pair[uint,double](udistance,distance)

    #Two dummy variables
    cdef double fabs_res = distance
    cdef size_t winner = udistance
    #Go over mutations in pop
    cdef unsigned twoN = 2*pop.popsize()
    for index in range(nmuts):
        #fwdpp allows extinct mutations to remain.  This is done
        #for efficiency reasons (fwdpp will recycle this spot in memory
        #.  We need to skip them here. We also skip any fixations.
        #(There will be no fixations in a typical "pop-gen" simulation,
        #but we have this check in here for completeness. During a typical
        #simualtion, fixations get moved from pop.mutations to pop.fixations,
        #and their fixation times get recorded in pop.fixation_times)  
        if pop.mcounts[index] > 0 and pop.mcounts[index] < twoN and pop.mutations[index].neutral is True:
            fabs_res = fabs(pop.mutations[index].pos - pos)
            if fabs_res < distance:
                #Mutation is closer to desired pos than previous
                if max(pop.mcounts[index],count) - min(pop.mcounts[index],count) < udistance:
                    #Mutation is closer to desired frequency than previos
                    winner = index
                    distance = fabs_res
                    udistance = max(pop.mcounts[index],count) - min(pop.mcounts[index],count)
    
    #We found it!
    pop.mutations[winner].s=new_s
    pop.mutations[winner].h=new_h
    #This next call is critical!!
    #Changing s,h is not sufficient to tell fwdpp
    #to use the mutation in fitness calculations.
    #We must call this fwdpp function to move the mutation
    #into the "selected variants" container of each gamete.
    change_neutral[singlepop_t](deref(pop),winner)
    return pair[uint,double](pop.mcounts[winner],pop.mutations[winner].pos)

def change_selection_coeff(PopType p,double pos, unsigned count, double new_s, double new_h):
    """
    This is our handy new Python function!

    It will find the NEUTRAL mutation in p whose position and number of occurrences is closest to pos and 
    count and change its s,h to new_s,new_h.

    If p is invariant, nothing will be done.
    """
    if isinstance(p,Spop):
        return change_s((<Spop>p).pop.get(),pos,count,new_s,new_h)
    else:
        raise RuntimeError("PopType not supported")
