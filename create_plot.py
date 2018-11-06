
import BNSVY2015
import CKS2018
import tiltedCHSH

import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import NullFormatter  # useful for `logit` scale





alpha=0


# plot the results from Bancal et al.
also_if_trivial = False
BNSVY2015_lower_bounds = BNSVY2015.lower_bounds(alpha=alpha,
                                                also_if_trivial=also_if_trivial)
BNSVY2015_violations = BNSVY2015.violations(alpha=alpha, also_if_trivial=also_if_trivial)
plt.plot(BNSVY2015_violations, BNSVY2015_lower_bounds, 'g.')

# plot the CKS18 results
CKS2018_violations = CKS2018.violations(alpha=alpha, also_if_trivial=False)
CKS2018_lower_bounds = \
        [ CKS2018.s(alpha=alpha) * beta + CKS2018.mu(alpha=alpha) for beta in CKS2018_violations]
plt.plot(CKS2018_violations, CKS2018_lower_bounds, 'r')

# plot trivial lower bound
if len(CKS2018_violations) > len(BNSVY2015_violations):
    violations = CKS2018_violations
else:
    violations = BNSVY2015_violations

trivial_lowbound = tiltedCHSH.trivial_lower_bound(alpha)
plt.plot(violations, np.array([trivial_lowbound for __ in violations]),'--',color='0.75')

# plot trivial upper bound
(a, b) = tiltedCHSH.trivial_upper_bound(alpha)
plt.plot(violations, [a * violation + b for violation in violations],':',color='0.75')
#
# plot layout
plt.title('hoi',fontsize=14)
#plt.xlim([classval(alpha),quantumval(alpha)+0.02])
#plt.ylim([ystart,1.03])
plt.xlabel(r'$\beta$',fontsize=14)

plt.show()

