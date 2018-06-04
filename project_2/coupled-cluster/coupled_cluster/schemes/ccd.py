import numpy as np

from .cc import CoupledCluster
from .cc_interface import amplitude_scaling_two_body

class CoupledClusterDoubles(CoupledCluster):

    def _initialize(self):
        o, v = self.o, self.v
        n, m, l = self.n, self.m, self.l

        self.f = self.h + np.einsum("piqi -> pq", self.u[:, o, :, o])

        self._compute_initial_guess()

        # Used for storing the previous value of self.t
        self._t = np.zeros(self.t.shape)

        # Used for intermediate terms
        self.term = np.zeros(self.t.shape)


        self.off_diag_f_bc = np.zeros((m, m))
        self.off_diag_f_kj = np.zeros((n, n))

        for b in range(m):
            for c in range(m):
                if b == c:
                    continue

                self.off_diag_f_bc[b, c] = self.f[b + n, c + n]

        for k in range(n):
            for j in range(n):
                if k == j:
                    continue

                self.off_diag_f_kj[k, j] = self.f[k, j]

    def _compute_initial_guess(self):
        o, v = self.o, self.v

        self.t = self.u[v, v, o, o].copy()
        amplitude_scaling_two_body(self.t, self.f, self.m, self.n)



    def _compute_ccd_energy(self):
        f, u, t, o, v = self.f, self.u, self.t, self.o, self.v

        e_ref = np.einsum("ii ->", f[o, o])
        e_ref -= 0.5 * np.einsum("ijij ->", u[o, o, o, o])
        e_ref += 0.25 * np.einsum("abij, abij ->", u[v, v, o, o], t)

        return e_ref

    def _compute_energy(self):
        return self._compute_ccd_energy()

    def _compute_amplitudes(self, theta):
        self._t.fill(0)

        self._compute_one_body_amplitude()
        self._compute_two_body_amplitude()

        amplitude_scaling_two_body(self._t, self.f, self.m, self.n)

        self.t = (1 - theta) * self._t + theta * self.t

    def _compute_one_body_amplitude(self):
        self.term = np.einsum(
                "bc, acij -> abij", self.off_diag_f_bc, self.t, out=self.term)
        self.term -= self.term.swapaxes(0, 1)
        self._t += self.term

        self.term = np.einsum(
                "kj, abik -> abij", self.off_diag_f_kj, self.t, out=self.term)
        self.term -= self.term.swapaxes(2, 3)
        self._t += self.term

    def _compute_two_body_amplitude(self):
        o, v = self.o, self.v

        self._t += self.u[v, v, o, o]

        self.term = 0.25 * np.einsum(
                "cdij, abmn, mncd -> abij", self.t, self.t, self.u[o, o, v, v],
                out=self.term)
        self._t += self.term

        self.term = 0.5 * np.einsum(
                "cdij, abcd -> abij", self.t, self.u[v, v, v, v], out=self.term)
        self._t += self.term

        self.term = 0.5 * np.einsum(
                "cdjm, abin, mncd -> abij", self.t, self.t, self.u[o, o, v, v],
                out=self.term)
        self.term -= self.term.swapaxes(2, 3)
        self._t += self.term

        self.term = - 0.5 * np.einsum(
                "acnm, bdij, nmcd -> abij", self.t, self.t, self.u[o, o, v, v],
                out=self.term)
        self.term -= self.term.swapaxes(0, 1)
        self._t += self.term

        self.term = np.einsum(
                "acim, bdjn, mncd -> abij", self.t, self.t, self.u[o, o, v, v],
                out=self.term)
        self.term -= self.term.swapaxes(2, 3)
        self._t += self.term

        self.term = np.einsum(
                "acim, bmjc -> abij", self.t, self.u[v, o, o, v], out=self.term)
        self.term -= self.term.swapaxes(0, 1)
        self.term -= self.term.swapaxes(2, 3)
        self._t += self.term

        self.term = 0.5 * np.einsum(
                "abmn, mnij -> abij", self.t, self.u[o, o, o, o], out=self.term)
        self._t += self.term
