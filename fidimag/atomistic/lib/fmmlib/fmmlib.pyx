from libcpp.vector cimport vector

cdef extern from "tree.hpp":
    cdef cppclass Particle:
        Particle(double x, double y, double z, double mux, double muy, double muz)
        double x, y, z, mux, muy, muz

    cdef cppclass Cell:
        Cell(double x, double y, double z, double r, unsigned int parent, unsigned int order,
            unsigned int level, unsigned int ncrit)
        unsigned int nleaf, nchild, level
        double x, y, z, r
        unsigned int parent
        vector[double] Mx
        vector[double] My
        vector[double] Mz
        vector[double] Lx
        vector[double] Ly
        vector[double] Lz
        vector[int] leaf

    cdef vector[Cell] build_tree(vector[Particle]& particles,
                             Cell &root,
                             unsigned int ncrit,
                             unsigned int order)


cdef vector[Particle] sim_to_particles(mesh, mu_s):
    """
    Convert the positions of the lattice spins to C++ Particle objects
    for the building of the trees.
    """
    pass



def build_tree(mesh, mus):
    pass
    # particles = 
    # cells = build_tree()




