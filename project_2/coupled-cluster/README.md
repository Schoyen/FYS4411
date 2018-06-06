# Coupled Cluster Doubles

This code uses a Python, Cython and C++ to run the coupled cluster doubles method on systems
of two-dimensional quantum dots.

## Installation

To install and compile the code run:
```bash
make install
```

This builds the project in the `/path/to/python/lib/python3.6/site-packages/` directory,
where /path/to/python/ depends on the default Python installation and `/python3.6/`
depends on the version of Python (Python 3.6 in this case).

To make sure that all the required libraries are installed and installing the project
you can run:
```bash
make installr
```

## Compilation

To compile the code in this directory run:
```bash
make build
```

## Testing

To check that the program is behaving the way it should you can run:
```bash
py.test
```
after the program has been installed/compiled. This will execute all the tests located in
[tests](coupled-cluster/tests).
