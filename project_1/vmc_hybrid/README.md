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
