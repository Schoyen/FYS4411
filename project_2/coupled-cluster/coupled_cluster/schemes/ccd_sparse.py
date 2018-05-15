import numpy as np
import sparse

from .cc import CoupledCluster

class CoupledClusterDoublesSparse(CoupledCluster):

    def _initialize(self, initial_guess):
        self.d = sparse.DOK((self.m, self.m, self.n, self.n))
        self._compute_d_matrix()

        self.h_sans_diag = sparse.DOK(self.h.shape)
        for p in range(len(self.h)):
            for q in range(len(self.h[p])):
                if p == q:
                    continue

                self.h_sans_diag[p, q] = self.h[p, q]

        self.h_sans_diag = self.h_sans_diag.to_coo()

        if initial_guess is not None:
            self.t = initial_guess
        else:
            self._compute_initial_guess()

        self.u_bc = sparse.DOK((self.m, self.m))
        self.u_kj = sparse.DOK((self.n, self.n))

        for b in range(self.m):
            for c in range(self.m):
                val = 0
                for k in range(self.n):
                    val += self.u[b + self.n, k, c + self.n, k]
                self.u_bc[b, c] = val

        self.u_bc = self.u_bc.to_coo()

        for k in range(self.n):
            for j in range(self.n):
                val = 0
                for l in range(self.n):
                    val += self.u[k, l, j, l]
                self.u_kj[k, j] = val

        self.u_kj = self.u_kj.to_coo()

    def _compute_d_matrix(self):
        h, n, m = self.h, self.n, self.m

        for a in range(m):
            for b in range(m):
                for i in range(n):
                    for j in range(n):
                        res = h[i, i] + h[j, j] \
                                - (h[n + a, n + a] + h[n + b, n + b])

                        if abs(res) < 1e-8:
                            continue

                        self.d[a, b, i, j] = 1.0/res

        self.d = self.d.to_coo()

    def _compute_initial_guess(self):
        u, d, o, v = self.u, self.d, self.o, self.v
        n = self.n

        self.t = u[v, v, o, o] * d

    def _compute_ccd_energy(self):
        h, u, t, o, v = self.h, self.u, self.t, self.o, self.v

        e_ref = np.einsum("ii->", h[o, o]) \
                + 0.5*np.einsum("ijij->", u[o, o, o, o])
        return e_ref + 0.25*np.einsum("abij, abij->", u[v, v, o, o], t)

    def _compute_energy(self):
        return self._compute_ccd_energy()

    def _compute_amplitudes(self):
        self._compute_intermediates()

        t_one_body = self._compute_one_body_amplitude()
        t_two_body = self._compute_two_body_amplitude()

        self.t = (t_one_body + t_two_body) * self.d

    def _compute_one_body_amplitude(self):
        h, t, o, v = self.h_sans_diag, self.t, self.o, self.v

        term = sparse.tensordot(
                h[o, o], t, axes=((0), (3))).transpose((1, 2, 0, 3))
        term_1 = term - term.transpose((0, 1, 3, 2))

        term = sparse.tensordot(
                h[v, v], t, axes=((1), (1))).transpose((0, 1, 3, 2))
        term_2 = term - term.transpose((1, 0, 2, 3))

        return term_1 + term_2

    def _compute_two_body_amplitude(self):
        h, u, t, o, v = self.h, self.u, self.t, self.o, self.v

        term = -sparse.tensordot(t, self.u_kj, axes=((3), (0)))
        term_1 = term - term.transpose((0, 1, 3, 2))

        term = sparse.tensordot(t, self.chi_lj, axes=((3), (0)))
        term_2 = term - term.transpose((0, 1, 3, 2))

        term_3 = sparse.tensordot(t, self.chi_klij, axes=((2, 3), (0, 1)))

        term = sparse.tensordot(t, self.chi_bc, axes=((1), (1)))
        term = term.transpose((0, 3, 2, 1)).transpose((0, 1, 3, 2))
        term_4 = term - term.transpose((1, 0, 2, 3))

        term = sparse.tensordot(t, self.chi_kbcj, axes=((1, 3), (2, 0)))
        term = term.transpose((0, 2, 1, 3))
        term = term - term.transpose((1, 0, 2, 3))
        term_5 = term - term.transpose((0, 1, 3, 2))

        term = 0.5*sparse.tensordot(t, u[v, v, v, v], axes=((0, 1), (2, 3)))
        term_6 = term.transpose((2, 3, 0, 1))

        return term_1 + term_2 + term_3 + term_4 + term_5 + term_6 \
                + u[v, v, o, o]

    def _compute_intermediates(self):
        h, u, t, o, v = self.h, self.u, self.t, self.o, self.v

        self.chi_klij = 0.25*sparse.tensordot(
                t, u[o, o, v, v], axes=((0, 1), (2, 3)))
        self.chi_klij = self.chi_klij.transpose((2, 3, 0, 1))
        self.chi_klij = self.chi_klij + 0.5*u[o, o, o, o]

        self.chi_bc = 0.5*sparse.tensordot(
                t, u[o, o, v, v], axes=((1, 2, 3), (2, 0, 1)))
        self.chi_bc = self.chi_bc + self.u_bc

        self.chi_lj = 0.5*sparse.tensordot(
                t, u[o, o, v, v], axes=((0, 1, 3), (2, 3, 0)))
        self.chi_lj = self.chi_lj.transpose()

        self.chi_kbcj = 0.5*sparse.tensordot(
                u[o, o, v, v], t, axes=((1, 3), (3, 1)))
        self.chi_kbcj = self.chi_kbcj.transpose((0, 2, 1, 3))
        self.chi_kbcj = self.chi_kbcj + u[o, v, v, o]
