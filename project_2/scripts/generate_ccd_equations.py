from sympy.physics.secondquant import (
        AntiSymmetricTensor, F, Fd, Commutator, wicks, evaluate_deltas,
        substitute_dummies, PermutationOperator, simplify_index_permutations, NO
)
from sympy import symbols, Dummy, Rational, factorial, latex

pretty_dummies = {
    'above': 'cdefgh',
    'below': 'mn',
    'general': 'pqrstu'
}

wicks_kwargs = {
    "simplify_dummies": True,
    "keep_only_fully_contracted": True,
    "simplify_kronecker_deltas": True
}

sub_kwargs = {
    "new_indices": True,
    "pretty_indices": pretty_dummies
}

i, j = symbols("i, j", below_fermi=True)
a, b = symbols("a, b", above_fermi=True)

def get_hamiltonian():
    p, q, r, s = symbols("p, q, r, s", cls=Dummy)
    f = AntiSymmetricTensor("f", (p,), (q,))
    u = AntiSymmetricTensor("u", (p, r), (q, s))

    f = f * NO(Fd(p) * F(q))
    u = u * NO(Fd(p) * Fd(q) * F(s) * F(r))

    return f, Rational(1, 4) * u

def get_doubles_cluster_operator():
    i, j = symbols("i, j", below_fermi=True, cls=Dummy)
    a, b = symbols("a, b", above_fermi=True, cls=Dummy)

    t = AntiSymmetricTensor("t", (a, b), (i, j))
    t = t * NO(Fd(a) * Fd(b) * F(j) * F(i))

    return [Rational(1, 4) * t]

def compute_hausdorff(f, cluster_func, num_terms=4):
    commutator = Commutator
    comm_term = f

    equation = comm_term

    for i in range(num_terms):
        t = sum(cluster_func())

        comm_term = wicks(commutator(comm_term, t))
        comm_term = substitute_dummies(evaluate_deltas(comm_term))

        equation += comm_term/factorial(i + 1)

    equation = equation.expand()
    equation = evaluate_deltas(equation)
    equation = substitute_dummies(
            equation, new_indices=True, pretty_indices=pretty_dummies)

    return equation


def get_energy_equation(equation_f, equation_u):
    eq = (equation_f + equation_u).expand()
    eq = evaluate_deltas(eq)
    eq = substitute_dummies(eq, **sub_kwargs)
    energy = wicks(eq, **wicks_kwargs)

    return energy

def get_one_body_equation(equation_f, equation_u):
    one_body_eq = wicks(
            NO(Fd(i) * Fd(j) * F(b) * F(a)) * equation_f, **wicks_kwargs)

    p = PermutationOperator
    one_body_eq = simplify_index_permutations(one_body_eq, [p(a, b), p(i, j)])
    one_body_eq = substitute_dummies(one_body_eq, **sub_kwargs)

    return one_body_eq

def get_two_body_equation(equation_f, equation_u):
    two_body_eq = wicks(
            NO(Fd(i) * Fd(j) * F(b) * F(a)) * equation_u, **wicks_kwargs)

    p = PermutationOperator
    two_body_eq = simplify_index_permutations(two_body_eq, [p(i, j), p(a, b)])
    two_body_eq = substitute_dummies(two_body_eq, **sub_kwargs)

    return two_body_eq

def get_ccd_equations():
    f, u = get_hamiltonian()
    equation_f = compute_hausdorff(f, get_doubles_cluster_operator)
    equation_u = compute_hausdorff(u, get_doubles_cluster_operator)

    energy = get_energy_equation(equation_f, equation_u)
    one_body = get_one_body_equation(equation_f, equation_u)
    two_body = get_two_body_equation(equation_f, equation_u)

    return energy, [one_body, two_body]

if __name__ == "__main__":
    energy, amplitudes = get_ccd_equations()
    print ("Energy rhs:")
    print (latex(energy))
    print ("One-body amplitude rhs:")
    print (latex(amplitudes[0]))
    print ("Two-body amplitude rhs:")
    print (latex(amplitudes[1]))
