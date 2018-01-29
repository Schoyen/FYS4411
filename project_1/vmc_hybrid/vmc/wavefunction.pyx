

cdef class Wavefunction:
    cdef parameters *m_parameters

    def __cinit__(self):
        self.m_parameters = get_variational_parameters()

    def __dealloc__(self):
        free_parameters_struct(self.m_parameters)
