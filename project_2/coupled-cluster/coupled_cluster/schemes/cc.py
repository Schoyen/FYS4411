import numpy as np

class CoupledCluster:

    def __init__(self, h, u, n, initial_guess=None):
        self.n = n
        self.l = len(h)
        self.m = self.l - self.n

        self.h = h
        self.u = u

        # Occupied (o) and virtual (v) array slices
        self.o = slice(0, self.n)
        self.v = slice(self.n, self.l)

        self._initialize(initial_guess=initial_guess)

    def __err(self):
        raise NotImplementedError("Use approximation subclass")

    def _initialize(self, initial_guess):
        self.__err()

    def _compute_energy(self):
        self.__err()

    def _compute_amplitudes(self, theta):
        self.__err()

    def compute_reference_energy(self):
        h, u, o, v = self.h, self.u, self.o, self.v
        e_ref = np.einsum("ii ->", h[o, o]) \
                + 0.5*np.einsum("ijij ->", u[o, o, o, o])

        return e_ref

    def compute_energy(self, max_iterations=100, tol=1e-4, theta=0.9):
        iterations = 0

        diff = 100
        energy = self._compute_energy()

        while diff > tol and iterations < max_iterations:
            self._compute_amplitudes(theta)
            energy_prev = energy
            energy = self._compute_energy()
            diff = abs(energy - energy_prev)
            iterations += 1

        return energy, iterations
