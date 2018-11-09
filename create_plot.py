
import BNSVY2015
import CKS2018
import tiltedCHSH

import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import NullFormatter  # useful for `logit` scale




def create_plot(ax, alpha, ystart, dotsize, num_decimals_beta_star):
    
    # plot the results from Bancal et al.
    also_if_trivial = False
    BNSVY2015_lower_bounds = BNSVY2015.lower_bounds(alpha=alpha,
                                                    also_if_trivial=also_if_trivial)
    BNSVY2015_violations = BNSVY2015.violations(alpha=alpha, also_if_trivial=also_if_trivial)
    ax.plot(BNSVY2015_violations, BNSVY2015_lower_bounds, 'g.', markersize=dotsize)
    
    # plot the CKS18 results
    CKS2018_violations = CKS2018.violations(alpha=alpha, also_if_trivial=False)
    CKS2018_lower_bounds = \
            [ CKS2018.s(alpha=alpha) * beta + CKS2018.mu(alpha=alpha) for beta in CKS2018_violations]
    ax.plot(CKS2018_violations, CKS2018_lower_bounds, 'r')
    
    # plot trivial lower bound
    violations = CKS2018.violations(alpha=alpha, also_if_trivial=True)
    
    trivial_lowbound = tiltedCHSH.trivial_lower_bound(alpha)
    ax.plot(violations, np.array([trivial_lowbound for __ in violations]),'--',color='0.75')
    
    # plot trivial upper bound
    (a, b) = tiltedCHSH.trivial_upper_bound(alpha)
    ax.plot(violations, [a * violation + b for violation in violations],':',color='0.75')


    # plot alpha
    cv = tiltedCHSH.classical_value(alpha)
    qv = tiltedCHSH.quantum_value(alpha)
    ax.text(cv + 0.14 * (qv - cv),
            1.03 - 0.4 * (1.03-ystart),
            r'$\alpha={}$'.format(alpha),
            fontsize=20)

    # plot beta-star
    ax.vlines(x=CKS2018.beta_star(alpha),
              ymin=ystart,
              ymax=tiltedCHSH.trivial_lower_bound(alpha),
              color='red',
              linewidths=1.0,
              zorder=2,
              linestyle='dashed') 
    ax.text(CKS2018.beta_star(alpha) + (qv - cv) * 0.02,
            (tiltedCHSH.trivial_lower_bound(alpha)-ystart)*0.35+ystart,
            r'$\beta^*\approx{}$'.format(np.round(CKS2018.beta_star(alpha),num_decimals_beta_star)),
            color='red',
            fontsize=18)


    ax.tick_params('both',labelsize=18)
    ax.set_ylim([ystart, 1.03])
    
    return ax




plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
params = {'text.usetex' : True,
        'font.size' : 20,
        'font.family' : 'lmodern',
        'text.latex.unicode' : True
        }
plt.rcParams.update(params)

axes = {}
fig, [axes[0], axes[0.5], axes[1]] = plt.subplots(3,1,figsize=(20,100))
fig.set_size_inches(6.5,17.5,forward=True)
for ax in list(axes.values()):
    ax.set_aspect('auto')


axes[1].set_xlabel('$\\beta_{\\alpha}$',fontsize=20)
plt.subplots_adjust(top=1.00, bottom=0.07, left=0.06, right=0.99)

ystarts = {0 : 0.4,
           0.5 : 0.6,
           1.0 : 0.75}
dotsizes = {0 : 4,
            0.5 : 7,
            1.0 : 7}

num_decimals_beta_star = {0 : 2,
                          0.5 : 2,
                          1.0 : 3}

for alpha in [0, 0.5, 1]:
    axes[alpha] = create_plot(ax=axes[alpha],
                              alpha=alpha,
                              ystart=ystarts[alpha],
                              dotsize=dotsizes[alpha],
                              num_decimals_beta_star=num_decimals_beta_star[alpha])



plt.show()
#plt.savefig("./output/CKS18_vs_BNSVY15.png")
