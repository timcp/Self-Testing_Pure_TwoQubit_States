import math
import numpy as np
import CKS2018
import matplotlib.pyplot as plt


def get_minimal_eigenvalue_of_hermitian(A):
    (eigenvalues, __ ) = np.linalg.eig(A)
    return np.real(min(eigenvalues))

def get_minimal_eigenvalue_of_hermitian_over_grid(A, grid, **additional_parameters):
    # set a start value
    a = grid[0][0]
    b = grid[0][1]
    current_minimum_eigenvalue = get_minimal_eigenvalue_of_hermitian(A(a=a, b=b, **additional_parameters))

    # compute the minimum eigenvalue
    for (a, b) in grid:
        m = get_minimal_eigenvalue_of_hermitian(A(a=a, b=b, **additional_parameters))
        if current_minimum_eigenvalue > m:
            current_minimum_eigenvalue = m
    return current_minimum_eigenvalue


if __name__ == "__main__":
    epsilon = 0.05 # grid coarseness
    
    alpha_start = 0.
    dalpha = 0.05
    alpha_stop = 2. - dalpha
    
    number_of_alphas = int(math.ceil( ( alpha_stop - alpha_start) / dalpha ) + 1.)
    alpha_list = np.linspace(start=alpha_start,
                             stop=alpha_stop,
                             num=number_of_alphas)
    
    
    
    # define the grid: discretization of [0, pi/4] x [0, pi/4]
    L = epsilon/np.sqrt(2.)
    number_of_grid_points = math.ceil(np.pi/(2. * L))
    
    x_points_list = np.linspace(0, np.pi/4., num=number_of_grid_points, endpoint=True)
    y_points_list = np.linspace(0, np.pi/2., num=number_of_grid_points, endpoint=True)
    
    grid = [(x, y) for x in x_points_list for y in y_points_list]
    
    
    
    
    print("Now computing the minimum eigenvalues for:")
    minevals = {}
    for alpha in alpha_list:
        print(r" +  \alpha = {} (will go up to but not including 2)".format(alpha))
        minevals[alpha] = get_minimal_eigenvalue_of_hermitian_over_grid(CKS2018.T,
                                                                        grid=grid,
                                                                        alpha=alpha)
    print(minevals)
    minimum = minevals[0]
    alpha_of_minimum = 0
    for alpha in alpha_list:
        m = minevals[alpha]
        if minimum > m:
            minimum = m
            alpha_of_minimum = alpha
    print("Minimum: {}".format(minimum))

    # store minimum eigenvalues
    description_string = \
"""
# This file contains values for alpha in [0, 2), the parameter of
# the family of tilted CHSH operators W_{alpha}, and the minimum 
# eigenvalue of the operator 
#   T_{alpha} = K_{alpha} - s_{alpha} * W_{alpha} - mu_{alpha} * Identity,
# minimized over a discretization of [0, pi/2] x [0, pi/4].
# Here, K_{alpha}, s_{alpha} and mu_{alpha} are defined as in the 
# article.
# 
# The first column contains values for alpha, the second the minimum
# eigenvalue of T_{alpha} over the grid.
#
# The minimum eigenvalue found over the entire grid is %s, which
# is achieved for alpha=%s.
#
"""%(minimum, alpha_of_minimum)

    data_save_location = "output/minimum_eigenvalues.csv"
    with open(data_save_location, "w+") as datafile:
        datafile.write(description_string)
        for alpha in alpha_list:
            datafile.write("{}, {}\n".format(alpha, minevals[alpha]))
    datafile.close()
    



    absminevals = [np.abs(mineval) for mineval in minevals]
    plt.yscale('log')
    plt.ylim([0, max(absminevals)*1.04])
    plt.scatter(alpha_list, absminevals)
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$|\min_{(a,b)}\quad\lambda_{\min} T_{\alpha}(a,b)|$")
    plot_save_location = "output/minimum_eigenvalues.pdf"
    plt.savefig(plot_save_location)
    print("A plot of the minimum eigenvalues can be found in {};\
           the data itself is stored in {}".format(plot_save_location, data_save_location))




