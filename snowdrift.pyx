from fwdpy.fwdpp cimport additive_diploid
from fwdpy.fitness cimport *
#Cython doesn't provide std::max?
cdef extern from "<algorithm>" namespace "std" nogil:
    T max[T](const T &, const T &)

    
cdef struct snowdrift_data:
    vector[double] phenotypes
    double b1,b2,c1,c2

cdef inline double snowdrift_fitness(const diploid_t & dip, const gcont_t & gametes, const mcont_t & mutations,
                              snowdrift_data & d) nogil:
    cdef double zself = d.phenotypes[<size_t>dip.e]
    cdef size_t i
    cdef double fitness = 0.0,zpair,a
    for i in range(d.phenotypes.size()):
        if i != <size_t>dip.e:
            zpair = zself + d.phenotypes[i]
            a = d.b1*zpair + d.b2*zpair - d.c1*zself - d.c2*zself*zself
            fitness += 1.0 + max[double](a,0.0)
    return fitness

cdef inline void snowdrift_update(singlepop_t * pop, snowdrift_data & d) nogil:
    if d.phenotypes.size() < pop.diploids.size():
        d.phenotypes.resize(pop.diploids.size())

    cdef size_t i=0
    cdef additive_diploid a
    for i in range(pop.diploids.size()):
        pop.diploids[i].e = <double>i
        d.phenotypes[i] = a(<diploid_t>pop.diploids[i],<gcont_t>pop.gametes,<mcont_t>pop.mutations,<double>2.0)

ctypedef singlepop_fitness_data[snowdrift_data] singlepop_fitness_snowdrift

cdef class SpopSnowdrift(SpopFitness):
    def __cinit__(self,double b1,double b2, double c1, double c2):
        cdef snowdrift_data d
        d.phenotypes = vector[double]()
        d.b1=b1
        d.b2=b2
        d.c1=c1
        d.c2=c2
        self.wfxn=singlepop_fitness_snowdrift(&snowdrift_fitness,
                                              &snowdrift_update)
