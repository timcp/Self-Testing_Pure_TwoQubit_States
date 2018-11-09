import tiltedCHSH
import math
import numpy as np
from tiltedCHSH import t, Id, X, Y, Z


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
    violations = tiltedCHSH.betas(alpha=alpha, numpoints=1000)
    if also_if_trivial:
        return violations
    else:
        bstar = beta_star(alpha)
        return [violation for violation in violations if violation >= bstar]


# defining the extraction channels

def b_star(alpha):
    return np.arcsin( np.sqrt( (4 - alpha*alpha) / 8. ) )

def Gamma(x,cutoffval):
    if x <= cutoffval:
        return X
    else:
        return Z

def conjugate(A, U):
    """
    Computes the matrix product :math:`UAU^{\dagger}`.
    """
    return np.matmul(np.matmul(U, A), U.getH())



def g(x):
    return (1.+np.sqrt(2.))*(np.cos(x) + np.sin(x) - 1.)


def effective_angle(alpha, x):
    bs = b_star(alpha)
    if x <= bs:
        return x * np.pi / (4.*bs)
    else:
        return np.pi/2. - np.pi/4. * (np.pi - 2*x)/(np.pi - 2*bs)


#def Lambda(angle, Gamma_operator, rhoAB, register=0):
#    g_value = g(angle)
#    if register == 0:
#        term = t(Gamma_operator,Id)
#    elif register == 1:
#        term = t(Id, Gamma_operator)
#    else:
#        raise Exception
#    return 0.5 * ((1. + g_value) * rhoAB +
#                  (1. - g_value) * term
#                 )
#
#def LambdaA(x, rhoAB):
#    return Lambda(angle=x,
#                  rhoAB=rhoAB,
#                  Gamma_operator=Gamma(x, np.pi/4.),
#                  register=0)
#
#
#def LambdaB(alpha,x,rhoAB):
#    return Lambda(angle=effective_angle(alpha=alpha, x=x),
#                  rhoAB=rhoAB,
#                  Gamma_operator=Gamma(x, b_star(alpha)),
#                  register=1)
           
def gB(alpha, x):
    return g(effective_angle(alpha=alpha, x=x))

def LambdaA(x,rhoAB):
    c0 = 0.5*(1.+g(x))
    c1 = 0.5*(1.-g(x))
    cutoffvalA = np.pi/4.
    return c0*rhoAB +\
            c1*conjugate(rhoAB,t(Gamma(x,cutoffvalA),Id))

def LambdaB(alpha,x,rhoAB):
    c0 = 0.5*(1.+gB(alpha,x))
    c1 = 0.5*(1.-gB(alpha,x))
    cutoffvalB = b_star(alpha)
    return c0*rhoAB +\
            c1*conjugate(rhoAB,t(Id,Gamma(x,cutoffvalB)))


# defining K and T
def K(alpha, a, b):
    return LambdaB(alpha, b, LambdaA(a, tiltedCHSH.Phi(alpha)))

def T(alpha, a, b):
    return K(alpha, a, b) - s(alpha) * tiltedCHSH.tilted_CHSH_operator(alpha, a, b) - mu(alpha) * t(Id, Id)


