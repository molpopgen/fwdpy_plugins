from fwdpy.fwdpp cimport additive_diploid
from fwdpy.fwdpy cimport Spop
from fwdpy.fitness cimport *
from libcpp.memory cimport unique_ptr
from libc.stdio cimport printf
cimport cython
#Cython doesn't provide std::max?
cdef extern from "<algorithm>" namespace "std" nogil:
    T max[T](const T &, const T &)

#fwdpy doesn't initialize the 'label' value of a diploid
#We do it here in case of situations where
#one also wants to set the population up with a specific mutation.
#Without this, we wouldn't a correct fitness calculation in the
#first generation
cdef inline void label_diploids_cpp(singlepop_t * pop):
    cdef size_t i=0
    for i in range(pop.diploids.size()):
        pop.diploids[i].label=i

#Data type to track stuff about trait values.
#The doubles are params of the snowdrift model
cdef cppclass snowdrift_data:
    vector[double] phenotypes
    double b1,b2,c1,c2
    site_dependent_fitness_wrapper a
    
#Function to return an initialized snowdrift_data object
cdef inline snowdrift_data make_snowdrift_data(double b1,double b2, double c1, double c2):
    cdef snowdrift_data d
    d.phenotypes=vector[double]()
    d.b1=b1
    d.b2=b2
    d.c1=c1
    d.c2=c2
    return d

#The fitness model:
cdef inline double snowdrift_fitness(const diploid_t & dip,
                                     const gcont_t & gametes,
                                     const mcont_t & mutations,
                                     snowdrift_data & d) nogil:
    cdef double zself = d.phenotypes[dip.label]
    cdef size_t i=0
    cdef double fitness = 0.0,zpair=0.0,a=0.0
    cdef size_t n = d.phenotypes.size()
    for i in range(n):
        if i != dip.label:
            zpair = zself + d.phenotypes[i]
            a = d.b1*zpair + d.b2*zpair*zpair - d.c1*zself - d.c2*zself*zself
            fitness += a
    fitness/=<double>(d.phenotypes.size()-1)
    return max[double](1+fitness,0.0)

#The updater function:
#Invidual fitness is additive over size, P = sum_i e_i,
#where e_i is i-th effect size
cdef inline void snowdrift_update(const singlepop_t * pop, snowdrift_data & d) nogil:
    if d.phenotypes.size() < pop.diploids.size():
        d.phenotypes.resize(pop.diploids.size())
    cdef size_t i=0
    for i in range(pop.diploids.size()):
        d.phenotypes[i] = d.a(<diploid_t>pop.diploids[i],
                              <gcont_t>pop.gametes,
                              <mcont_t>pop.mutations,
                              <genotype_fitness_updater>hom_additive_update_2,
                              <genotype_fitness_updater>het_additive_update,
                              <fitness_function_finalizer>return_trait_value,
                              <double>0.0)

#Nicer name for our stateful fitness type:
ctypedef singlepop_fitness_data[snowdrift_data] singlepop_fitness_snowdrift

cdef class SpopSnowdrift(SpopFitness):
    """
    This is the fitness model as a Cython extension class
    """
    def __cinit__(self,double b1,double b2, double c1, double c2):
        """
        Constructor takes the params for the payoff functions.
        """
        cdef snowdrift_data d = make_snowdrift_data(b1,b2,c1,c2)
        self.wfxn=<unique_ptr[singlepop_fitness]>unique_ptr[singlepop_fitness_snowdrift](new singlepop_fitness_snowdrift(snowdrift_fitness,
                                                                                                                         snowdrift_update,d))

def label_diploids(Spop pop):
    """
    Python fxn to label the diploids
    """
    label_diploids_cpp(pop.pop.get())
