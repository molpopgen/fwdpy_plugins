#Example of plugins implementing "stateless" fitness models
#We will implement the standard additive model several different
#ways, plus one haplotype-based model

#The additive model of fitness is 1+sum(g_i), where g_i is sh
#for Aa genotypes and 2s for aa genotypes.

#The additive model of a trait value would simply be sum(g_i),
#and therefore the "return fitness" function matters:

from fwdpy.fitness cimport *
from libcpp.memory cimport unique_ptr

cdef extern from "<algorithm>" namespace "std" nogil:
    T max[T](T,T)

#First, we'll implement the additive model
#by writing functions that process Aa and aa 
#genotypes at each site and then return a final fitness 
#value

#This function is how fitness is updated for heterozygous (Aa) 
#genotypes
cdef inline void additiveAa(double & w, const popgenmut & m) nogil:
    (&w)[0] += (m.h*m.s)

#This function is how fitness is updated for a homozygous genotype
cdef inline void additiveaa(double & w, const popgenmut & m) nogil:
    (&w)[0] += (2.0*m.s)

#For additive models, we can simply return 1+w.
#If we wanted additivity for a quantitative trait,
#then this function would return w, to 
#center the value at an average trait value of 0.
cdef inline double return_fitness(double w) nogil:
    return max[double](1.0+w,0.0)

#Now, we construct our custom fitness class,
#which can be used as a fitness model in a simulation
cdef class CustomAdditive1(SpopFitness):
    def __cinit__(self):
        self.wfxn = unique_ptr[singlepop_fitness](new singlepop_fitness(additiveAa,     #Pass in the "Aa updating fxn"
									additiveaa,     #The "aa updating fxn"
								    	return_fitness, #The fxn to return the final fitness value
								  	0.0))           #The initial value for fitness

#See fwdpy/fitness.pxd and fwdpy/fitness.pyx for how additive & multiplicative fitness models are implemented in fwdpy

#Now, let's implement it as a more complex function representing the full signature of a fitness function:

#import fwdpp's function object for calculating additve effects over sites:
from fwdpy.fwdpp cimport additive_diploid

cdef inline double additive_fitness_fxn1(const diploid_t & dip, const gcont_t & gametes, const mcont_t & mutations) nogil:
    #Create an additive_diploid object:
    cdef additive_diploid a
    return a(dip,gametes,mutations,2.0)

cdef inline double additive_fitness_fxn2(const diploid_t & dip, const gcont_t & gametes, const mcont_t & mutations) nogil:
    #Use the wrapper around fwdpp's site_dependent fitness
    #defined in fwdpy/fitness.pxd:
    cdef site_dependent_fitness_wrapper wrapper
    #call it like a function, passing our cdef functions from above:
    return wrapper(dip,gametes,mutations,additiveaa,additiveAa,return_fitness,0.0) #Note that additiveaa is passed in b4 additiveAa.  Yes, this is inconsistent w.r.to above, but it matches the fwdpp interface for site_depdendent_fitness.

#Now we can define two new additive fitness models:

cdef class CustomAdditive2(SpopFitness):
    def __cinit__(self):
        self.wfxn = unique_ptr[singlepop_fitness](new singlepop_fitness(additive_fitness_fxn1))

cdef class CustomAdditive3(SpopFitness):
    def __cinit__(self):
        self.wfxn = unique_ptr[singlepop_fitness](new singlepop_fitness(additive_fitness_fxn2))

#Finally, let's define additive fitness as the sum of the sums of selection coefficients on each
#hapltype.  Note that this explicitly assumes h=1!!!!

#We will use functions from fwdpy/fitness.pxd in order to get the job done:

cdef double add_hap_effects(double a,double b) nogil:
    return max[double](1.0+a+b,0.0)

cdef class CustomAdditive4(SpopFitness):
    def __cinit__(self):
        self.wfxn=unique_ptr[singlepop_fitness](new singlepop_fitness(sum_haplotype_effects,add_hap_effects))
