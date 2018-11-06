import math


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
    return math.sqrt( 8 + 2 * alpha)


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


def betas(alpha):
    """auxillary function"""
    qval = quantum_value(alpha=alpha)
    start = 2.0
    delta = 0.01
    violations = []
    beta = start
    while beta < qval:
        violations.append(beta)
        beta += delta
    return violations



