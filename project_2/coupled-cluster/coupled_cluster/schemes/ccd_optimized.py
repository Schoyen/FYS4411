import numpy as np
import time

from .ccd import CoupledClusterDoubles
from .cc_interface import amplitude_scaling_two_body
from .helper_ccd import (
        compute_chi_abcd, compute_chi_bmjc, compute_chi_ad, compute_chi_nj,
        compute_chi_abcd_contraction, compute_f_bc_t_contraction,
        compute_f_kj_t_contraction, compute_chi_bmjc_contraction,
        compute_t_u_contraction, compute_chi_nj_contraction,
        compute_chi_ad_contraction
)

class CoupledClusterDoublesOptimized(CoupledClusterDoubles):

    def _initialize(self, **kwargs):
        super(CoupledClusterDoublesOptimized, self)._initialize(**kwargs)

        o, v = self.o, self.v

        # Allocate memory for the intermediates
        self.chi_abcd = np.zeros((self.m, self.m, self.m, self.m))
        self.chi_ad = np.zeros((self.m, self.m))
        self.chi_bmjc = np.zeros((self.m, self.n, self.n, self.m))
        self.chi_bjcm = np.zeros((self.m, self.n, self.m, self.n))
        self.chi_nj = np.zeros((self.n, self.n))

        self._t_ijab = np.zeros((self.n, self.n, self.m, self.m))
        self._t_iabj = np.zeros((self.n, self.m, self.m, self.n))
        self._t_aibj = np.zeros((self.m, self.n, self.m, self.n))
        self._u_hpph = np.zeros((self.n, self.m, self.m, self.n))
        self._u_phph = np.zeros((self.m, self.n, self.m, self.n))
        self._u_phph_2 = np.zeros((self.m, self.n, self.m, self.n))

        self._t_ijab[:] = self.t.transpose(2, 3, 0, 1)
        self._t_iabj[:] = self.t.transpose(2, 0, 1, 3)
        self._t_aibj[:] = self.t.transpose(0, 2, 1, 3)
        self._u_hpph[:] = self.u[o, o, v, v].transpose(0, 2, 3, 1)
        self._u_phph[:] = self.u[o, o, v, v].transpose(2, 0, 3, 1)
        self._u_phph_2[:] = self.u[v, o, o, v].transpose(0, 2, 3, 1)

    def _compute_amplitudes(self, theta):
        self._t.fill(0)

        #if self.parallel:
        #    self._compute_intermediates_parallel()
        #    self._compute_one_body_amplitude_parallel()
        #    self._compute_two_body_amplitude_parallel()
        #else:
        #    self._compute_intermediates()
        #    self._compute_one_body_amplitude()
        #    self._compute_two_body_amplitude()

        print ("Starting parallel timing...")
        self._compute_intermediates_parallel()
        self._compute_one_body_amplitude_parallel()
        self._compute_two_body_amplitude_parallel()

        print ("\nStarting serial timing...")
        self._t.fill(0)
        self._compute_intermediates()
        self._compute_one_body_amplitude()
        self._compute_two_body_amplitude()

        print ("\n")

        amplitude_scaling_two_body(self._t, self.f, self.m, self.n)
        self.t = np.add((1 - theta) * self._t, theta * self.t, out=self.t)

        self._t_ijab[:] = self.t.transpose(2, 3, 0, 1)
        self._t_iabj[:] = self.t.transpose(2, 0, 1, 3)
        self._t_aibj[:] = self.t.transpose(0, 2, 1, 3)

    @profile
    def _compute_one_body_amplitude_parallel(self):
        t0 = time.time()
        compute_f_bc_t_contraction(
                self.term, self.off_diag_f_bc, self.t, self.n, self.m)
        t1 = time.time()
        print ("parallel f_bc_t contraction: {0} sec".format(t1 - t0))
        self._t += self.term

        t0 = time.time()
        compute_f_kj_t_contraction(
                self.term, self.off_diag_f_kj, self.t, self.n, self.m)
        t1 = time.time()
        print ("parallel f_kj_t contraction: {0} sec".format(t1 - t0))
        self._t += self.term

    @profile
    def _compute_two_body_amplitude_parallel(self):
        o, v = self.o, self.v

        self._t += self.u[v, v, o, o]

        t0 = time.time()
        compute_chi_abcd_contraction(
                self.term, self._t_ijab, self.chi_abcd, self.n, self.m)
        t1 = time.time()
        print ("parallel chi_abcd contraction: {0} sec".format(t1 - t0))
        self._t += self.term

        t0 = time.time()
        compute_chi_nj_contraction(
                self.term, self.t, self.chi_nj, self.n, self.m)
        t1 = time.time()
        print ("parallel chi_nj contraction: {0} sec".format(t1 - t0))
        self._t += self.term

        t0 = time.time()
        compute_chi_ad_contraction(
                self.term, self.t, self.chi_ad, self.n, self.m)
        t1 = time.time()
        print ("parallel chi_ad contraction: {0} sec".format(t1 - t0))
        self._t += self.term

        t0 = time.time()
        compute_chi_bmjc_contraction(
                self.term, self.t, self.chi_bmjc, self.n, self.m)
        t1 = time.time()
        print ("parallel chi_bmjc contraction: {0} sec".format(t1 - t0))
        self._t += self.term

        t0 = time.time()
        compute_t_u_contraction(self.term, self.t, self.u, self.n, self.m)
        t1 = time.time()
        print ("parallel t_u contraction: {0} sec".format(t1 - t0))
        self._t += self.term

    @profile
    def _compute_intermediates_parallel(self):
        o, v = self.o, self.v

        t0 = time.time()
        compute_chi_abcd(
                self.chi_abcd, self.t, self.u[o, o, v, v],
                self.u[v, v, v, v], self.n, self.m)
        t1 = time.time()
        print ("parallel chi_abcd: {0} sec".format(t1 - t0))

        t0 = time.time()
        compute_chi_ad(self.chi_ad, self.t, self.u, self.n, self.m)
        t1 = time.time()
        print ("parallel chi_ad: {0} sec".format(t1 - t0))

        t0 = time.time()
        compute_chi_bmjc(
                self.chi_bmjc, self._t_iabj, self._u_hpph,
                self.u[v, o, o, v], self.n, self.m)
        t1 = time.time()
        print ("parallel chi_bmjc: {0} sec".format(t1 - t0))

        t0 = time.time()
        compute_chi_bjcm(
                self.chi_bjcm, self._t_aibj, self._u_phph,
                self._u_phph_2, self.n, self.m)
        t1 = time.time()
        print ("parallel chi_bmjc: {0} sec".format(t1 - t0))

        t0 = time.time()
        compute_chi_nj(self.chi_nj, self.t, self.u, self.n, self.m)
        t1 = time.time()
        print ("parallel chi_nj: {0} sec".format(t1 - t0))

    @profile
    def _compute_one_body_amplitude(self):
        t0 = time.time()
        self.term = np.einsum(
                "bc, acij -> abij", self.off_diag_f_bc, self.t,
                out=self.term, optimize="optimal")
        self.term -= self.term.swapaxes(0, 1)
        t1 = time.time()
        print ("f_bc_t contraction: {0} sec".format(t1 - t0))
        self._t += self.term

        t0 = time.time()
        self.term = np.einsum(
                "kj, abik -> abij", self.off_diag_f_kj, self.t,
                out=self.term, optimize="optimal")
        self.term -= self.term.swapaxes(2, 3)
        t1 = time.time()
        print ("f_kj_t contraction: {0} sec".format(t1 - t0))
        self._t += self.term

    @profile
    def _compute_two_body_amplitude(self):
        o, v = self.o, self.v

        self._t += self.u[v, v, o, o]

        t0 = time.time()
        self.term = np.einsum(
                "cdij, abcd -> abij", self.t, self.chi_abcd,
                out=self.term, optimize="optimal")
        t1 = time.time()
        print ("chi_abcd contraction: {0} sec".format(t1 - t0))
        self._t += self.term

        t0 = time.time()
        self.term = np.einsum(
                "abin, nj -> abij", self.t, self.chi_nj,
                out=self.term, optimize="optimal")
        self.term -= self.term.swapaxes(2, 3)
        t1 = time.time()
        print ("chi_nj contraction: {0} sec".format(t1 - t0))
        self._t += self.term

        t0 = time.time()
        self.term = - np.einsum(
                "bdij, ad -> abij", self.t, self.chi_ad,
                out=self.term, optimize="optimal")
        self.term -= self.term.swapaxes(0, 1)
        t1 = time.time()
        print ("chi_ad contraction: {0} sec".format(t1 - t0))
        self._t += self.term

        t0 = time.time()
        self.term = np.einsum(
                "acim, bmjc -> abij", self.t, self.chi_bmjc,
                out=self.term, optimize="optimal")
        self.term -= self.term.swapaxes(0, 1)
        self.term -= self.term.swapaxes(2, 3)
        t1 = time.time()
        print ("chi_bmjc contraction: {0} sec".format(t1 - t0))
        self._t += self.term

        t0 = time.time()
        self.term = 0.5 * np.einsum(
                "abmn, mnij -> abij", self.t, self.u[o, o, o, o],
                out=self.term, optimize="optimal")
        t1 = time.time()
        print ("t_u contraction: {0} sec".format(t1 - t0))
        self._t += self.term

    @profile
    def _compute_intermediates(self):
        o, v = self.o, self.v

        t0 = time.time()
        self.chi_abcd = 0.25 * np.einsum(
                "abmn, mncd -> abcd", self.t, self.u[o, o, v, v],
                out=self.chi_abcd, optimize="optimal")
        self.chi_abcd += 0.5 * self.u[v, v, v, v]
        t1 = time.time()
        print ("chi_abcd: {0} sec".format(t1 - t0))

        t0 = time.time()
        self.chi_ad = 0.5 * np.einsum(
                "acnm, nmcd -> ad", self.t, self.u[o, o, v, v],
                out=self.chi_ad, optimize="optimal")
        t1 = time.time()
        print ("chi_ad: {0} sec".format(t1 - t0))

        t0 = time.time()
        self.chi_bmjc = 0.5 * np.einsum(
                "bdjn, mncd -> bmjc", self.t, self.u[o, o, v, v],
                out=self.chi_bmjc, optimize="optimal")
        self.chi_bmjc += self.u[v, o, o, v]
        t1 = time.time()
        print ("chi_bmjc: {0} sec".format(t1 - t0))

        t0 = time.time()
        self.chi_nj = 0.5 * np.einsum(
                "cdjm, mncd -> nj", self.t, self.u[o, o, v, v],
                out=self.chi_nj, optimize="optimal")
        t1 = time.time()
        print ("chi_nj: {0} sec".format(t1 - t0))
