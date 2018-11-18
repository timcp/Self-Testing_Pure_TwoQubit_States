#####################################################
#                                                   #
# This file contains a number of Python methods     #
# that define the quantities used in the proof of   #
# Conjecture 1, such as s_{\alpha}, \mu_{\alpha},   #
# the extraction channels and the operators K       #
# and T as given on page 15.                        #
#                                                   #
# This file is part of the supporting material      #
# belonging to:                                     #
# "Robust self-testing of two-qubit states"         #
# Tim Coopmans, Jędrzej Kaniewski and Christian     #
# Schaffner (2018)                                  #
# arXiv:                                            #
#                                                   #
#####################################################


import tiltedCHSH
import numpy as np
from tiltedCHSH import t, Id, X, Y, Z


#####################################################
# Definitions of quantities that occur in the main  #
# text of the article.                              #
#####################################################

def s(alpha):
    """
    Returns the parameter :math:`s` as given in Conjecture 1
    of the article.
    """
    sqrt1 = np.sqrt( 4 - alpha * alpha) / tiltedCHSH.quantum_value(alpha)
    sqrt2 = np.sqrt( 2 * alpha * alpha) / tiltedCHSH.quantum_value(alpha)
    return ( 1 - (1 + sqrt1 + sqrt2 )/4. ) / (tiltedCHSH.quantum_value(alpha)-tiltedCHSH.classical_value(alpha))

def mu(alpha):
    """
    Returns the parameter :math:`\mu` as given in Conjecture 1
    of the article.
    """
    return 1 - s(alpha) * tiltedCHSH.quantum_value(alpha)

def beta_star(alpha):
    """
    Returns the infimum of all violations, denoted as :math:`\beta^*`,
    for which the lower bound to the extractability, as given in Conjecture 1
    of the article, exceeds the trivial lower bound.
    
    To be precise, :math:`\beta` is the solution to the equation
    ``` s_{\alpha} \cdot \beta + \mu_{\alpha} = \lambda_{\alpha}^2```
    where :math:`\lambda_{\alpha}^2` is the square of the largest
    Schmidt coefficient of the target state and thus the trivial lower 
    bound to the extractability. See also the last equation on page 3
    of the article.
    """
    return ( tiltedCHSH.trivial_lower_bound(alpha) - mu(alpha) ) / s(alpha)

def b_star(alpha):
    """
    The optimal angle on Bob's side, as given just below eq. (11)
    in the article.
    """
    return np.arcsin( np.sqrt( (4 - alpha*alpha) / 8. ) )


#####################################################
# Definitions of quantities that occur in           #
# Appendix B of the article.                        #
#####################################################


def Gamma(x):
    if x <= np.pi/4.:
        return X
    else:
        return Z



def g(x):
    return ( 1. + np.sqrt(2.) ) * ( np.cos(x) + np.sin(x) - 1. )

def effective_angle(alpha, x):
    """
    Denoted in Appendix B as :math:`h_{\alpha}`.
    """
    bs = b_star(alpha)
    if x <= bs:
        return x * np.pi / (4.*bs)
    else:
        return np.pi/2. - np.pi/4. * (np.pi - 2*x)/(np.pi - 2*bs)

def conjugate(A, U):
    """
    Auxillary function.
    Computes the matrix product :math:`UAU^{\dagger}`.
    """
    return np.matmul(np.matmul(U, A), U.getH())

def _Lambda(angle, rhoAB, register=0):
    """
    Auxillary function.

    Returns the output of Alice's (Bob's) extraction
    channel on the bipartite state `rhoAB` in case
    the parameter `register` is set to 0 (1).

    To be precise: this method returns the output of
    the channel
    .. math ::
            \rho \mapsto (\Lambda_A(x) \otimes I)(\rho)
    
    on input state `rhoAB`, where :math:`\Lambda_A(x)(\rho)`
    is given on page 15 of the article and :math:`I` is
    the identity channel. In case the parameter
    `register` is set to 1, then this method returns
    the output of the channel
    .. math ::
            \rho \mapsto (I \otimes \Lambda_B(x))(\rho)

    on input state `rhoAB`.
    """

    # Define on which of the two registers
    # the channel should act
    if register == 0:
        U = t(Gamma(angle), Id)
    elif register == 1:
        U = t(Id, Gamma(angle))
    else:
        raise Exception

    g_value = g(angle)
    return 0.5 * ((1. + g_value) * rhoAB +
                  (1. - g_value) * conjugate(rhoAB, U)
                 )

def LambdaA(x, rhoAB):
    """
    The output of Alice's extraction channel
    on the bipartite input state `rhoAB`.
    """
    return _Lambda(angle=x,
                   rhoAB=rhoAB,
                   register=0)


def LambdaB(alpha, x, rhoAB):
    """
    The output of Bob's extraction channel
    on the bipartite input state `rhoAB`.
    """
    return _Lambda(angle=effective_angle(alpha=alpha, x=x),
                   rhoAB=rhoAB,
                   register=1)

def K(alpha, a, b):
    """
    The target state after both Alice and Bob applied
    their extraction channels. Also stated on page 15
    of the article.

    Parameters
    ----------
    alpha : float in the interval [0, 2)
        Parameter that determines the target state.
    a : float
        Alice's angle.
    b : float
        Bob's angle.
    """
    return LambdaB(alpha, b, LambdaA(a, tiltedCHSH.target_state(alpha)))

def T(alpha, a, b):
    """
    The operator as given on page 15 of the article.
    """
    return K(alpha, a, b)\
           - s(alpha) * tiltedCHSH.tilted_CHSH_operator(alpha, a, b)\
           - mu(alpha) * t(Id, Id)
