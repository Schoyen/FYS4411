import numpy as np

from .ccd import CoupledClusterDoubles

class CoupledClusterDoublesOptimized(CoupledClusterDoubles):

    def _initialize(self):
        super(CoupledClusterDoublesOptimized, self)._initialize()

    def _compute_amplitudes(self, theta):
        # TODO: Use DIIS

        self._t.fill(0)
