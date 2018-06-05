import numpy as np

from .ccd import CoupledClusterDoubles
from .cc_interface import amplitude_scaling_two_body
from .helper_ccd import (
        compute_chi_abcd, compute_chi_bmjc, compute_chi_abcd_contraction,
        compute_chi_bmjc_contraction, compute_t_u_contraction,
        compute_chi_ad, compute_chi_nj
)

class CoupledClusterDoublesOptimized(CoupledClusterDoubles):

    def _initialize(self, **kwargs):
        super(CoupledClusterDoublesOptimized, self)._initialize(**kwargs)

        # Allocate memory for the intermediates
        self.chi_abcd = np.zeros((self.m, self.m, self.m, self.m))
        self.chi_ad = np.zeros((self.m, self.m))
        self.chi_bmjc = np.zeros((self.m, self.n, self.n, self.m))
        self.chi_nj = np.zeros((self.n, self.n))

    def _compute_amplitudes(self, theta):
        self._t.fill(0)

        if self.parallel:
            self._compute_intermediates_parallel()
            self._compute_one_body_amplitude()
            self._compute_two_body_amplitude_parallel()
        else:
            self._compute_intermediates()
            self._compute_one_body_amplitude()
            self._compute_two_body_amplitude()

        amplitude_scaling_two_body(self._t, self.f, self.m, self.n)
        self.t = np.add((1 - theta) * self._t, theta * self.t, out=self.t)

    def _compute_two_body_amplitude_parallel(self):
        o, v = self.o, self.v

        self._t += self.u[v, v, o, o]

        compute_chi_abcd_contraction(
                self.term, self.t, self.chi_abcd, self.n, self.m)
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

        compute_chi_bmjc_contraction(
                self.term, self.t, self.chi_bmjc, self.n, self.m)
        self._t += self.term

        compute_t_u_contraction(self.term, self.t, self.u, self.n, self.m)

        self._t += self.term

    def _compute_intermediates_parallel(self):
        o, v = self.o, self.v

        compute_chi_abcd(self.chi_abcd, self.t, self.u, self.n, self.m)
        compute_chi_ad(self.chi_ad, self.t, self.u, self.n, self.m)
        compute_chi_bmjc(self.chi_bmjc, self.t, self.u, self.n, self.m)
        compute_chi_nj(self.chi_nj, self.t, self.u, self.n, self.m)

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
