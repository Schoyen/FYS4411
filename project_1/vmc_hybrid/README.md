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

## Implementation tests
We currently reproduce the results, i.e., the plots, from [these slides](https://www.acsu.buffalo.edu/~phygons/cp2/topic5/topic5-lec1.pdf).

## To do

All figures __MUST__ be on the xkcd-format.

### Initial testing
#### Tables
Include different variances.

1. Numerical and analytical results for one dimension and one particle. CZHECH
2. Numerical and analytical results for `num_samples = [10, 100, 500]` and `num_dimensions = [2, 3]` for `alpha = 0.5`.

#### Figures

1. Figures corresponding to tables in appendices. Compare with exact answer.
2. Variance plots. Compare variance ratios.

### Importance sampling
Its important to discuss how much faster (in terms of number of samples) we reach an equilibrium state.

#### Tables
1. Compare results (energy and variance) as a function of `step_length` for a fixed `alpha`.

#### Figures
1. Compare variance decrease as a function of cycles for importance sampling and brute force.
2. Try to plot variance as a function of `step_length` for different alphas.

### Repulsive interaction

#### Tables
1. Compute statistics for a set of alphas.
2. Compute statistics for "true alpha" found from steepest descent.
  - Benchmark with references.

#### Figures
1. Plot energy and variance for table 1.
2. Plot of alphas (show convergence for several alphas) found from steepest descent.

### Onebody densities
Cool figures with and without Jastrow.
