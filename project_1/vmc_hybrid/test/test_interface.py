from vmc.interface import PySimpleGaussian, PyHarmonicOscillator
import numpy as np
import pytest

@pytest.fixture
def wavefunction_params():
    spread = 2.0
    num_particles = 100
    num_parameters = 1
    num_dimensions = 3

    return spread, num_particles, num_parameters, num_dimensions

@pytest.fixture()
def one_simple_gauss():
    spread = 1.0
    num_particles = 1
    num_parameters = 1
    num_dimensions = 1

    return spread, num_particles, num_parameters, num_dimensions

@pytest.fixture()
def two_simple_gauss():
    spread = 1.0
    num_particles = 1
    num_parameters = 1
    num_dimensions = 2

    return spread, num_particles, num_parameters, num_dimensions

@pytest.fixture()
def three_simple_gauss():
    spread = 1.0
    num_particles = 1
    num_parameters = 1
    num_dimensions = 3

    return spread, num_particles, num_parameters, num_dimensions

@pytest.fixture
def hamiltonian_params():
    mass = 1
    omega = 1

    return mass, omega

def test_wavefunction(wavefunction_params):
    spread, num_particles, num_parameters, num_dimensions = wavefunction_params

    parameters = np.array([0.5])

    wavefunction = PySimpleGaussian(
            num_particles, num_dimensions, num_parameters, spread=spread)

    particles = wavefunction.get_particles()
    assert ((particles >= -spread) & (particles <= spread)).all()

    wavefunction.set_parameters(parameters)
    assert np.allclose(parameters, wavefunction.get_parameters())


@pytest.mark.run(after="test_wavefunction")
def test_one_dimensional_ho(hamiltonian_params, one_simple_gauss):
    mass, omega = hamiltonian_params
    spread, num_particles, num_parameters, num_dimensions = one_simple_gauss

    ho = PyHarmonicOscillator(mass, omega)
    wavefunction = PySimpleGaussian(
            num_particles, num_dimensions, num_parameters, spread=spread)

    parameters = np.array([0.5])
    wavefunction.set_parameters(parameters)

    local_energy = ho.compute_local_energy(wavefunction)
    expected_value = 0.5

    assert np.allclose(local_energy, expected_value)

@pytest.mark.run(after="test_one_dimensional_ho")
def test_two_dimensional_ho(hamiltonian_params, two_simple_gauss):
    mass, omega = hamiltonian_params
    spread, num_particles, num_parameters, num_dimensions = two_simple_gauss

    ho = PyHarmonicOscillator(mass, omega)
    wavefunction = PySimpleGaussian(
            num_particles, num_dimensions, num_parameters, spread=spread)

    parameters = np.array([0.5])
    wavefunction.set_parameters(parameters)

    local_energy = ho.compute_local_energy(wavefunction)
    expected_value = 1.0

    assert np.allclose(local_energy, expected_value)

@pytest.mark.run(after="test_two_dimensional_ho")
def test_three_dimensional_ho(hamiltonian_params, three_simple_gauss):
    mass, omega = hamiltonian_params
    spread, num_particles, num_parameters, num_dimensions = three_simple_gauss

    ho = PyHarmonicOscillator(mass, omega)
    wavefunction = PySimpleGaussian(
            num_particles, num_dimensions, num_parameters, spread=spread)

    parameters = np.array([0.5])
    wavefunction.set_parameters(parameters)

    local_energy = ho.compute_local_energy(wavefunction)
    expected_value = 1.5

    assert np.allclose(local_energy, expected_value)