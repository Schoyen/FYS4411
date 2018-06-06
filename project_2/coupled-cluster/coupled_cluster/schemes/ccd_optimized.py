import numpy as np

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
        n, m = self.n, self.m

        # Allocate memory for the intermediates
        self.chi_abcd = np.zeros((self.m, self.m, self.m, self.m))
        self.chi_ad = np.zeros((self.m, self.m))
        # We use different intermediates for the parallel situation and the
        # serial implementation as the latter requires the intermediate on the
        # form of a matrix.
        if self.parallel:
            self.chi_bmjc = np.zeros((self.m, self.n, self.n, self.m))
        else:
            self.chi_cmbj = np.zeros((m * n, m * n))
        self.chi_nj = np.zeros((self.n, self.n))

        # Allocate memory for blocks
        if self.parallel:
            self._t_ijab = self.t.transpose(2, 3, 0, 1).copy()
            self._t_iabj = self.t.transpose(2, 0, 1, 3).copy()
            self._u_hpph = self.u[o, o, v, v].transpose(0, 2, 3, 1).copy()
        else:
            self.u_hhhh = self.u[o, o, o, o]\
                    .reshape(n * n, n * n)
            self.u_hhpp = self.u[o, o, v, v]\
                    .reshape(n * n, m * m)
            self.u_hhpp_phh_p = self.u[o, o, v, v]\
                    .transpose(2, 0, 1, 3)\
                    .reshape(m * n * n, m)
            self.u_hhpp_phph = self.u[o, o, v, v]\
                    .transpose(3, 1, 2, 0)\
                    .reshape(m * n, m * n)
            self.u_phhp_phph = self.u[v, o, o, v]\
                    .transpose(0, 2, 3, 1)\
                    .reshape(m * n, m * n)
            self.u_hhpp_h_pph = self.u[o, o, v, v]\
                    .transpose(1, 2, 3, 0)\
                    .reshape(n, m * m * n)

            self.t_phph = self.t\
                    .transpose(0, 2, 1, 3)\
                    .reshape(m * n, m * n)
            self.t_pph_h = self.t\
                    .transpose(0, 1, 3, 2)\
                    .reshape(m * m * n, n)
            self.t_p_phh = self.t\
                    .transpose(1, 0, 2, 3)\
                    .reshape(m, m * n * n)

    def _compute_amplitudes(self, theta):
        self._t.fill(0)

        if self.parallel:
            self._compute_intermediates_parallel()
            self._compute_one_body_amplitude_parallel()
            self._compute_two_body_amplitude_parallel()
        else:
            self._compute_intermediates()
            self._compute_one_body_amplitude()
            self._compute_two_body_amplitude()

        amplitude_scaling_two_body(self._t, self.f, self.m, self.n)
        self.t = np.add((1 - theta) * self._t, theta * self.t, out=self.t)

        # Update block forms of the amplitudes
        if self.parallel:
            self._t_ijab[:] = self.t.transpose(2, 3, 0, 1)
            self._t_iabj[:] = self.t.transpose(2, 0, 1, 3)
        else:
            self.t_phph[:] = self.t\
                    .transpose(0, 2, 1, 3)\
                    .reshape(self.m * self.n, self.m * self.n)
            self.t_pph_h[:] = self.t\
                    .transpose(0, 1, 3, 2)\
                    .reshape(self.m * self.m * self.n, self.n)
            self.t_p_phh[:] = self.t\
                    .transpose(1, 0, 2, 3)\
                    .reshape(self.m, self.m * self.n * self.n)

    def _compute_one_body_amplitude_parallel(self):
        compute_f_bc_t_contraction(
                self.term, self.off_diag_f_bc, self.t, self.n, self.m)
        self._t += self.term

        compute_f_kj_t_contraction(
                self.term, self.off_diag_f_kj, self.t, self.n, self.m)
        self._t += self.term

    def _compute_two_body_amplitude_parallel(self):
        o, v = self.o, self.v

        self._t += self.u[v, v, o, o]

        compute_chi_abcd_contraction(
                self.term, self._t_ijab, self.chi_abcd, self.n, self.m)
        self._t += self.term

        compute_chi_nj_contraction(
                self.term, self.t, self.chi_nj, self.n, self.m)
        self._t += self.term

        compute_chi_ad_contraction(
                self.term, self.t, self.chi_ad, self.n, self.m)
        self._t += self.term

        compute_chi_bmjc_contraction(
                self.term, self.t, self.chi_bmjc, self.n, self.m)
        self._t += self.term

        compute_t_u_contraction(self.term, self.t, self.u, self.n, self.m)
        self._t += self.term

    def _compute_intermediates_parallel(self):
        o, v = self.o, self.v

        compute_chi_abcd(
                self.chi_abcd, self.t, self.u[o, o, v, v],
                self.u[v, v, v, v], self.n, self.m)

        compute_chi_ad(self.chi_ad, self.t, self.u, self.n, self.m)

        compute_chi_bmjc(
                self.chi_bmjc, self._t_iabj, self._u_hpph,
                self.u[v, o, o, v], self.n, self.m)

        compute_chi_nj(self.chi_nj, self.t, self.u, self.n, self.m)

    def _compute_one_body_amplitude(self):
        m, n = self.m, self.n

        np.matmul(
                self.off_diag_f_bc,
                self.t_p_phh,
                out=self.term.reshape(m, m * n * n)
        )
        self.term[:] = self.term.transpose(1, 0, 2, 3)
        self.term -= self.term.swapaxes(0, 1)
        self._t += self.term

        np.matmul(
                self.t.reshape(m * m * n, n),
                self.off_diag_f_kj,
                out=self.term.reshape(m * m * n, n)
        )
        self.term *= -1
        self._t += self.term

    def _compute_two_body_amplitude(self):
        o, v = self.o, self.v
        n, m = self.n, self.m

        self._t += self.u[v, v, o, o]

        np.matmul(
                self.chi_abcd.reshape(m * m, m * m),
                self.t.reshape(m * m, n * n),
                out=self.term.reshape(m * m, n * n)
        )
        self._t += self.term

        np.matmul(
                self.t.reshape(m * m * n, n),
                self.chi_nj,
                out=self.term.reshape(m * m * n, n)
        )
        self.term -= self.term.swapaxes(2, 3)
        self._t += self.term

        np.matmul(
                self.chi_ad,
                self.t_p_phh,
                out=self.term.reshape(m, m * n * n)
        )
        self.term *= -1
        self.term -= self.term.swapaxes(0, 1)
        self._t += self.term

        self.term[:] = np.matmul(
                self.t_phph,
                self.chi_cmbj
        ).reshape(m, n, m, n).transpose(0, 2, 1, 3)
        self.term -= self.term.swapaxes(0, 1)
        self.term -= self.term.swapaxes(2, 3)
        self._t += self.term


        np.matmul(
                self.t.reshape(m * m, n * n),
                self.u_hhhh,
                out=self.term.reshape(m * m, n * n)
        )
        self.term *= 0.5
        self._t += self.term

    def _compute_intermediates(self):
        o, v = self.o, self.v
        n, m = self.n, self.m

        np.matmul(
                self.t.reshape(m * m, n * n),
                self.u_hhpp,
                out=self.chi_abcd.reshape(m * m, m * m)
        )
        self.chi_abcd *= 0.25
        self.chi_abcd += 0.5 * self.u[v, v, v, v]

        np.matmul(
                self.t.reshape(m, m * n * n),
                self.u_hhpp_phh_p,
                out=self.chi_ad
        )
        self.chi_ad *= 0.5

        np.matmul(
                self.t_phph,
                self.u_hhpp_phph,
                out=self.chi_cmbj
        )
        self.chi_cmbj *= 0.5
        self.chi_cmbj += self.u_phhp_phph
        self.chi_cmbj[:] = self.chi_cmbj.transpose(1, 0)

        np.matmul(
                self.u_hhpp_h_pph,
                self.t_pph_h,
                out=self.chi_nj
        )
        self.chi_nj *= 0.5
