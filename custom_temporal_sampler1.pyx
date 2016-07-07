#This is an example of a custom temporal sampler.

# Support for custom samplers is limited in fwdpy, but we are working on it.
# The cimport statements below from fwdpy.fwdpy show the minimial stuff required
# to implement a custom sampler for a single-deme simulation. uint is not needed in general,
# but is only here for this example.

from fwdpy.fwdpy cimport TemporalSampler,sampler_base,custom_sampler,singlepop_t,uint
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.memory cimport unique_ptr

#How do custom temporal sampler work?  They rely on defining an instance of the C++ template class
#custom_sampler.  The template type will be the return value from the sampler.  That return value will
#be a container of the recorded observations over time.

#This specific example records the number of occurences and selected coefficient of every selected mutation
#over time.  In c++, that is a pair<unsigned,double>, where the unsigned value is the count, and the double is
#the selection coefficient.  For each generation, there is a vector of such pairs.  Let's define a typedef
#for these per-generation observations:

ctypedef vector[pair[uint,double]] observation_t

#Ok, we need to associate the observation_t with a generation.  Let's use another pair:
                                                            
ctypedef pair[uint,observation_t] generation_data_t

#Finally, a vector of such things will be sampled over time:

ctypedef vector[generation_data_t] final_t

#We call that "final_t" because that will be the final type collected by a sampler for a single replicate simulation.
#The get() function of our sampler will return a vector[final_t], representing the sampled data for  all populations
#that we simulate.

#Now, we can define our sampler type:
ctypedef custom_sampler[final_t] SelectedSFSSampler_t

#This is the function that will actually do the sampling.
#For a custom_sampler[FOO], the "signature" of the function must be:
# void(*)(const singlepop_t *, const unsigned, FOO &):
cdef void update_selected_sfs(const singlepop_t * pop, const unsigned generation,
                              final_t & rv) nogil:
    #This function is reponsible for filling up rv with generation_data_t objects
    cdef observation_t data
    cdef size_t i=0
    cdef size_t n=pop.mutations.size()
    cdef unsigned twoN=2*pop.popsize()
    for i in range(n):
        if pop.mcounts[i]>0 and pop.mcounts[i] < twoN and pop.mutations[i].neutral is False:
            data.push_back(pair[uint,double](pop.mcounts[i],pop.mutations[i].s))
    rv.push_back(pair[uint,vector[pair[uint,double]]](generation,data))
    
cdef class SelectedSFSSampler(TemporalSampler):
    """
    This is our Cython extension class that will be used in Python.
    """
    def __cinit__(self,unsigned n):
        """
        Constructor

        :param n: A length.  Must be equal to the length of the SpopVec that is is being applied to
        """
        for i in range(n):
            #This is a Cython annoyance alert!
            #In C++11, if you have base class A and derived class B, a unique_ptr<B> is IMPLICITLY CONVERTIBLE
            #to unique_ptr<A>.  Cython doesn't recognize that, but this typecast (e.g., <A>B) suffices.
            #KEY POINT: we pass a function pointer to update_selected_sfs to the constructor of our custom sampler type
            self.vec.push_back(<unique_ptr[sampler_base]>unique_ptr[SelectedSFSSampler_t](new SelectedSFSSampler_t(&update_selected_sfs)))
    def get(self):
        """
        This collects the data from each element, creates a vector[sample_t], and returns it.

        This "just works" because final_t is a vector<pair<unsigned,vector<pair<unsigned,double> > >.

        Let's break that down:

        pair<unsigned,double> will be implicitly converted to a Python tuple.
        vector<pair<unsigned,double> > will be implicitly converted to a list of tuples.
        etc.

        Thus, this function returns a list.  That list contains the data for each generation in the form of a tuple.
        The first member of each tuple is the generation.  The second member is a list of tuples representing mutation count
        and selection coefficient over time.

        In our example Python script, we're going to clean up these data and coerce them into a pandas.DataFrame.

        If we weren't so lazy, we may convert this whole thing to a list of C++ maps, which are implicitly convertible to dicts,
        which are directly coercable to Pandas DataFrames.
        """
        cdef vector[final_t] rv
        cdef size_t i=0;
        cdef n=self.vec.size()
        for i in range(n):
            rv.push_back((<SelectedSFSSampler_t*>self.vec[i].get()).final())
        return rv

    
