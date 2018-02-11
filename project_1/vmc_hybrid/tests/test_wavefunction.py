from vmc.interface import PySimpleGaussian
import numpy as np


def test_wavefunction():
    spread = 2.0
    num_particles = 100
    num_parameters = 1
    num_dimensions = 3

    parameters = np.array([0.5])

    wavefunction = PySimpleGaussian(
            num_particles, num_dimensions, num_parameters, spread=spread)

    particles = wavefunction.get_particles()
    assert ((particles >= -spread) & (particles <= spread)).all()

    wavefunction.set_parameters(parameters)
    assert np.allclose(parameters, wavefunction.get_parameters())
