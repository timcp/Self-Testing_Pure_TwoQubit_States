import math
import numpy as np

# some matrices
Id = np.matrix([[1., 0.], [0., 1.]])
X = np.matrix([[0., 1.], [1., 0.]])
Y = np.matrix([[0., -1.j], [1.j, 0.]])
Z = np.matrix([[1., 0.],[0., -1.]])
H = (X + Z) / np.sqrt(2)
Hminus = (X - Z) / np.sqrt(2)

def t(A, B):
    return np.kron(A, B)

def Phi(alpha):
    """
    The target state, a partially-entangled
    two-qubit state.

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

def obs(r, x): #observable
    return np.cos(x) * X + (1. - 2. * r) * np.sin(x) * Z

def CHSH_operator(a, b):
    terms = [ (1 - 2. * k * m) * t(obs(k, a), obs(m, b))\
                for (k, m) in [(0, 0),
                               (0, 1),
                               (1, 0),
                               (1, 1)]]
    return sum(terms)


def tilted_CHSH_operator(alpha, a, b):
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
    a = (1-trivial_lower_bound(alpha)) / (quantum_value(alpha)-classical_value(alpha)) 
    b = 1 - a * quantum_value(alpha) 
    return (a,b)


def betas(alpha, numpoints=None):
    """auxillary function"""
    qval = quantum_value(alpha=alpha)
    beta = classical_value(alpha=alpha)
    violations = []
    if numpoints is None:
        delta = 0.01
    else:
        delta = (qval - beta)/numpoints
    while beta < qval:
        violations.append(beta)
        beta += delta
    return violations



