#####################################################
#                                                   #
# This file contains a number of Python methods     #
# that define the family of tilted CHSH operators   #
# and related quantities, such as its classical     #
# value, its quantum value and the state that it    #
# self-tests.                                       #
#                                                   #
# This file is part of the supporting material      #
# belonging to:                                     #
# "Robust self-testing of two-qubit states"         #
# Tim Coopmans, Jędrek Kaniewski and Christian      #
# Schaffner (2018)                                  #
# arXiv:                                            #
#                                                   #
#####################################################

import math
import numpy as np

# The Pauli and Hadamard matrices
Id  = np.matrix([[1.,   0.], [0., 1.]])
X   = np.matrix([[0.,   1.], [1., 0.]])
Y   = np.matrix([[0., -1.j], [1.j, 0.]])
Z   = np.matrix([[1.,   0.], [0., -1.]])
H   = (X + Z) / np.sqrt(2)
Hminus = (X - Z) / np.sqrt(2)



def t(A, B):
    """
    Abbreviated notation for the Kronecker product
    (tensor product) of two matrices.
    """
    return np.kron(A, B)

def target_state(alpha):
    """
    The target state that the tilted-CHSH-operator
    self-tests. The target state is, up to local unitaries,
    equal to :math:`\cos(\theta)|00\rangle + \sin(\theta)|11\rangle`.

    The value :math:`\alpha \in [0, 2)` "sweeps out"
    all possible pure partially-entangled two-qubit states,
    up to local unitary transformation.
    """
    return \
            0.25*\
            (\
            t(Id,Id) + \
            + np.sqrt( (2*alpha*alpha)/(4 + alpha*alpha) ) * (t(H,Id) + t(Id,X))\
            + t(H, X)
            + np.sqrt( (4 - alpha*alpha) / (4 + alpha*alpha) ) * (t(Y,Y) + t(Hminus, Z))\
            )

def obs(r, angle):
    """
    The binary observables as given in eq. (8) and (9) in the 
    article.
    """
    return np.cos(angle) * X + (1. - 2. * r) * np.sin(angle) * Z

def CHSH_operator(a, b):
    """
    The CHSH operator, defined as function of two angles,
    one on Alice's side and one on Bob's.
    """
    terms = [ (1 - 2. * k * m) * t(obs(k, a), obs(m, b))\
                for (k, m) in [(0, 0),
                               (0, 1),
                               (1, 0),
                               (1, 1)]]
    return sum(terms)

def tilted_CHSH_operator(alpha, a, b):
    """
    The tilted-CHSH-operator, parametrized by :math:`alpha`
    and the two angles, one for each party.
    Also given in eq. (10) of the article.
    """
    return alpha * t(obs(0, a), Id) + CHSH_operator(a, b)

def classical_value(alpha):
    """
    Returns the classical value of the tilted CHSH inequality
    with parameter :math:`\alpha \in [0, 2)`.
    """
    return 2 + alpha

def quantum_value(alpha):
    """
    Returns the quantum value of the tilted CHSH inequality
    with parameter :math:`\alpha \in [0, 2)`.
    """
    return math.sqrt( 8 + 2 * alpha * alpha)

def trivial_lower_bound(alpha):
    """
    Returns the trivial lower bound to the extractability of the
    state :math:`\cos(\theta_{\alpha}) |00\rangle + \sin(\theta_{\alpha} |11\rangle`
    where :math:`\theta_{\alpha}` is given in the paper.
    The parameter :math:`\alpha \in [0, 2)`.
    """
    return 0.5 * ( math.sqrt( ( 2 * alpha * alpha / (4 + alpha * alpha) ) ) + 1)

def trivial_upper_bound(alpha):
    """
    Computes the trivial upper bound to the extractability of the
    state :math:`\cos(\theta_{\alpha}) |00\rangle + \sin(\theta_{\alpha} |11\rangle`
    where :math:`\theta_{\alpha}` is given in the paper.
    The parameter :math:`\alpha \in [0, 2)`.

    Returns
    -------
    Tuple (float, float)
        (Slope, point where line crosses the vertical axis)
    """
    a = (1 - trivial_lower_bound(alpha)) / (quantum_value(alpha) - classical_value(alpha))
    b = 1 - a * quantum_value(alpha) 
    return (a,b)

def possible_violations(alpha, number_of_points=None):
    """
    Auxillary function. Discretizes the interval
    :math:`[c_{\alpha}, q_{\alpha}]`, where
    :math:`c_{\alpha}` is the classical value
    and :math:`q_{\alpha}` the quantum value of
    the tilted-CHSH-operator.

    Returns
    -------
    List of equally-spaced values, ranging
    from classical value to quantum value 
    of the tilted-CHSH-operator.
    """
    qval = quantum_value(alpha=alpha)
    beta = classical_value(alpha=alpha)
    violations = []
    if number_of_points is None:
        delta = 0.01
    else:
        delta = (qval - beta) / number_of_points
    while beta < qval:
        violations.append(beta)
        beta += delta
    return violations
