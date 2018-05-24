import numpy as np
import scipy.linalg

def calculate_energy(F, W, D):
    J = np.einsum("sr, qrps -> qp", D, W)
    K = np.einsum("sr, qrsp -> qp", D, W)

    energy = np.einsum("pq, qp ->", F, D)
    energy -= 0.5*np.einsum("pq, qp ->", J, D)
    energy += 0.25*np.einsum("pq, qp ->", K, D)

    return energy

def build_density_matrix(U, num_occupied):
    """Function building a density matrix.

    Build U_occ consisting of K x N (rows x columns), where K is the number of
    basis functions in total and N is the number of occupied basis functions.

    The prodcut np.dot(U_occ, U_occ.conj().T) then builds a K x K matrix, which
    is the density matrix for a given spin direction.

    Args:
        U: A two-dimensional NumPy ndarray found by solving the generalized
            eigenvalue problem for the Roothan-Hall and/or the Pople-Nesbet
            equations.
        num_occupied: The number of occupied states in the current
            spin-direction.

    Returns:
        np.ndarray: The density matrix as a two-dimensional NumPy ndarray.
    """

    U_occ = U[:, 0:num_occupied]
    return np.dot(U_occ, U_occ.conj().T)

def build_rhf_fock_matrix(H, W, D):
    """Function building the Fock matrix.

    The elements of the Fock matrix is given by:

        F_{pq} = H_{pq} + J(D)_{pq} - 0.5*K(D)_{pq},

    where we build J and K from the electron repulsion integrals (W) (using
    Mulliken notation) and the density matrix D.

        J(D)_{pq} = (pq|rs)*D_{sr} = W[p][q][r][s]*D[s][r],
        K(D)_{pq} = (ps|rq)*D_{sr} = W[p][s][r][q]*D[s][r].

    The NumPy-function einsum does this in a form which resemble the
    mathematical-product.

    Args:
        H: A two-dimensional NumPy ndarray with the one-body integrals as
            elements.
        W: A four-dimensional NumPy ndarray with the two-body electron
            repulsion integrals as elements.
        D: The density matrix of the system as a two-dimensional NumPy ndarray.

    Returns:
        np.ndarray: The RHF-Fock matrix.
    """
    J = np.einsum("sr, qrps -> qp", D, W)
    K = np.einsum("sr, qrsp -> qp", D, W)

    return H + J - 0.5*K

def scf_rhf(H, W, S, num_occupied, theta=0.01, tol=1e-8, max_iterations=1000):
    """The SCF-scheme for solving the RHF Roothan-Hall equations.

    Args:
        H: A two-dimensional NumPy ndarray with the one-body integrals as
            elements.
        W: A four-dimensional NumPy ndarray with the two-body electron
            repulsion integrals as elements.
        S: A two-dimensional NumPy ndarray with the overlap between the atomic
            orbital basis as elements.
        num_occupied: The number of occupied atomic orbitals (counting both spin
            directions).
        theta: A parameter to toggle mixing from the previous density matrix
            into the current one.
        tol: Value defining the convergence criteria for the SCF-iterations.
        iteration_warning: A boolean-value used when the user wants the program
            to get interrupted in case of divergence.

    Returns:
        float: The RHF-energy sans the nuclear repulsion energy.
    """

    # Make sure that theta is a valid number
    assert (0 <= theta <= 1), ("Theta should be in [0, 1]")

    # Initialize a counter to check if the iteration diverges or gets stuck
    counter = 0

    # Solve the initial generalized eigenvalue problem using the one-body
    # integrals as an approximation to the Fock matrix and the overlap matrix
    energy, U = scipy.linalg.eigh(H, S)

    # Get the density matrix from the occupied states in the U-matrix
    D = 2*build_density_matrix(U, num_occupied)
    # Build the Fock matrix from this system
    F = build_rhf_fock_matrix(H, W, D)

    # Set initial energy diff
    diff = 100

    # Loop until the solution converges, i.e., until the change in the
    # eigen-energies are lower than the threshold tol
    while diff > tol and counter < max_iterations:
        # Set the previous energy
        energy_prev = energy

        # Calculate the new density matrix adding in a theta-fraction of the
        # previous density matrix (if the need arises)
        D = (1 - theta) * 2 * build_density_matrix(U, num_occupied) \
            + theta * D

        # Build the updated Fock matrix
        F = build_rhf_fock_matrix(H, W, D)
        # Solve the new generalized eigenvalue problem
        energy, U = scipy.linalg.eigh(F, S)

        # Compute new difference
        diff = np.max(np.abs(energy - energy_prev))

        # Increment the counter
        counter += 1

    # Return the conversion matrix U
    return U, calculate_energy(F, W, D)
