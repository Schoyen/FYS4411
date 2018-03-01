# Variational Monte Carlo

This code uses a hybrid solution with Python, Cython and C++ to run Quantum
Variational Monte Carlo simulations.

## Compilation

To compile the code use:

```bash
python3 setup.py build_ext
```

This builds a shared library which, by default, gets placed in the `vmc`
directory.

## Installation
To do a global installation run

```bash
pip install . --upgrade
```

This builds the project in the `/path/to/python/lib/python3.6/site-packages/`
directory, where `/path/to/python/` depends on the default Python installation
and `/python3.6/` depends on the version of Python (Python 3.6 in this case).

### Local installation
An alternative to this can be to define your own local build directory which can
later be purged. For example:

```bash
mkdir -p $HOME/build/bin
mkdir -p $HOME/build/lib/python3.6/site-packages
```

We then need to set the `PATH` and `PYTHONPATH` environment variables.

```bash
export PATH=$PATH:$HOME/build/bin
export PYTHONPATH=$PYTHONPATH:$HOME/build/lib/python3.6/site-packages
```

Installation can then be performed by executing

```bash
pip install . --upgrade --install-options="--prefix=$HOME/build"
```

## To be implemented

Below follows a list of things that needs to be implemented/tested.

1. Implement the Metropolis sampling function in `src/metropolis_sampling.c`.
2. Add bosonic hard sphere wavefunction, variational parameters, local energy
   and ratio functions in `src/bosonic_hard_sphere.c`.
3. Test the implementation for the simple one-dimensional system described in
   Thijssen.
4. Use `valgrind` to locate memory leaks.

## Implementation tests
We currently reproduce the results, i.e., the plots, from [these slides](https://www.acsu.buffalo.edu/~phygons/cp2/topic5/topic5-lec1.pdf).
