import tiltedCHSH
import math

def s(alpha):
    """
    Returns the parameter :math:`s` in the lower bound 
    :math:`s_{\alpha} \cdot \beta + \mu_{\alpha}`
    as function of the violation :math:`\beta`, as proven in [CKS18].
    """
    sqrt1 = math.sqrt( 4 - alpha * alpha) / tiltedCHSH.quantum_value(alpha)
    sqrt2 = math.sqrt( 2 * alpha * alpha) / tiltedCHSH.quantum_value(alpha)
    return ( 1 - (1 + sqrt1 + sqrt2 )/4. ) / (tiltedCHSH.quantum_value(alpha)-tiltedCHSH.classical_value(alpha))

def mu(alpha):
    """
    Returns the parameter :math:`\mu` in the lower bound 
    :math:`\mu_{\alpha} \cdot \beta + \mu_{\alpha}`
    as function of the violation :math:`\beta`, as proven in [CKS18].
    """
    return 1 - s(alpha) * tiltedCHSH.quantum_value(alpha)

def beta_star(alpha):
    """
    Returns the infimum of all violations, denoted as :math:`\beta^*`,
    for which the lower bound to the extractability as given in [CKS18]
    exceeds the trivial lower bound.
    
    To be precise, :math:`\beta` is the solution to the equation
    ``` s_{\alpha} \cdot \beta + \mu_{\alpha} = {\ell}_{\alpha}```
    where :math:`\ell_{\alpha}` is the trivial lower bound to the extractability.
    and the line :math:`s_{\alpha} \cdot \beta + \mu_{\alpha}` as function of 
    the violation :math:`\beta` is a lower bound to the extractability.
    """
    return ( tiltedCHSH.trivial_lower_bound(alpha) - mu(alpha) ) / s(alpha)


def violations(alpha, also_if_trivial=True):
    violations = tiltedCHSH.betas(alpha=alpha, numpoints=100)
    if also_if_trivial:
        return violations
    else:
        bstar = beta_star(alpha)
        return [violation for violation in violations if violation >= bstar]


