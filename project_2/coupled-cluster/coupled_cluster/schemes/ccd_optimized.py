import numpy as np

from .ccd import CoupledClusterDoubles
from .cc_interface import amplitude_scaling_two_body

class CoupledClusterDoublesOptimized(CoupledClusterDoubles):

    def _initialize(self):
        super(CoupledClusterDoublesOptimized, self)._initialize()

        # Allocate memory for the intermediates
        self.chi_abcd = np.zeros((self.m, self.m, self.m, self.m))
        self.chi_ad = np.zeros((self.m, self.m))
        self.chi_bmjc = np.zeros((self.m, self.n, self.n, self.n))
        self.chi_nj = np.zeros((self.n, self.n))

    def _compute_amplitudes(self, theta):
        # TODO: Use DIIS

        self._compute_intermediates()
        super(CoupledClusterDoublesOptimized, self)._compute_amplitudes(theta)


        self._compute_one_body_amplitude()
        self._compute_two_body_amplitude()

    def _compute_one_body_amplitude(self):
        pass

    def _compute_two_body_amplitude(self):
        pass

    def _compute_intermediates(self):
        o, v = self.o, self.v

        self.chi_abcd = 0.25 * np.einsum(
                "abmn, mncd -> abcd", self.t, self.u[o, o, v, v],
                out=self.chi_abcd, optimize="optimal")
        self.chi_abcd += 0.5 * self.u[v, v, v, v]

        self.chi_ad = 0.5 * np.einsum(
                "acnm, nmcd -> ad", self.t, self.u[o, o, v, v],
                out=self.chi_ad, optimize="optimal")

        self.chi_bmjc = 0.5 * np.einsum(
                "bdjn, mncd -> bmjc", self.t, self.u[o, o, v, v],
                out=self.chi_bmjc, optimize="optimal")
        self.chi_bmjc += self.u[v, o, o, v]

        self.chi_nj = 0.5 * np.einsum(
                "cdjm, mncd -> nj", self.t, self.u[o, o, v, v],
                out=self.chi_nj, optimize="optimal")
