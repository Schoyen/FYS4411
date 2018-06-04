import numpy as np

from .ccd import CoupledClusterDoubles
from .cc_interface import amplitude_scaling_two_body

class CoupledClusterDoublesOptimized(CoupledClusterDoubles):

    def _initialize(self, max_diis_size=100, **kwargs):
        super(CoupledClusterDoublesOptimized, self)._initialize(**kwargs)

        # Allocate memory for the intermediates
        self.chi_abcd = np.zeros((self.m, self.m, self.m, self.m))
        self.chi_ad = np.zeros((self.m, self.m))
        self.chi_bmjc = np.zeros((self.m, self.n, self.n, self.m))
        self.chi_nj = np.zeros((self.n, self.n))

    def _compute_amplitudes(self, theta):
        # TODO: Use DIIS

        self._compute_intermediates()
        super(CoupledClusterDoublesOptimized, self)._compute_amplitudes(theta)

    def _compute_two_body_amplitude(self):
        o, v = self.o, self.v

        self._t += self.u[v, v, o, o]

        self.term = np.einsum(
                "cdij, abcd -> abij", self.t, self.chi_abcd,
                out=self.term, optimize="optimal")
        self._t += self.term

        self.term = np.einsum(
                "abin, nj -> abij", self.t, self.chi_nj,
                out=self.term, optimize="optimal")
        self.term -= self.term.swapaxes(2, 3)
        self._t += self.term

        self.term = - np.einsum(
                "bdij, ad -> abij", self.t, self.chi_ad,
                out=self.term, optimize="optimal")
        self.term -= self.term.swapaxes(0, 1)
        self._t += self.term

        self.term = np.einsum(
                "acim, bmjc -> abij", self.t, self.chi_bmjc,
                out=self.term, optimize="optimal")
        self.term -= self.term.swapaxes(0, 1)
        self.term -= self.term.swapaxes(2, 3)
        self._t += self.term

        self.term = 0.5 * np.einsum(
                "abmn, mnij -> abij", self.t, self.u[o, o, o, o],
                out=self.term, optimize="optimal")
        self._t += self.term

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
