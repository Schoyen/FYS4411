import numpy as np

from .cc import CoupledCluster
from .cc_interface import _amplitude_scaling_two_body

class CoupledClusterDoubles(CoupledCluster):

    def _initialize(self):
        self._compute_initial_guess()

    def _compute_initial_guess(self):
        o, v = self.o, self.v

        self.t = self.u[v, v, o, o]
        _amplitude_scaling_two_body(self.t, self.h, self.m, self.n)
