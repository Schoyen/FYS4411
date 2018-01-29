

cdef class MetropolisSampling:

    def __cinit__(self):
        m_wavefunction = Wavefunction()
        m_particles = Particles()
