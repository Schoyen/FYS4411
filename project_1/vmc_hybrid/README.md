# Variational Monte Carlo

This code uses a hybrid solution with Python, Cython and C to run Quantum
Variational Monte Carlo simulations.

## Compilation

To compile the code use:

```bash
python3 setup.py build_ext
```

This builds a shared library which, by default, gets placed in the `scripts`
directory. Scripts located in this directoy can import directly from this
library.

## To be implemented

Below follows a list of things that needs to be implemented/tested.

1. Implement the Metropolis sampling function in `src/metropolis_sampling.c`.
2. Add bosonic hard sphere wavefunction, variational parameters, local energy
   and ratio functions in `src/bosonic_hard_sphere.c`.
3. Test the implementation for the simple one-dimensional system described in
   Thijssen.
4. Use `valgrind` to locate memory leaks.
