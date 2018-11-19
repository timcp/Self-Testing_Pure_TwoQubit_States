#####################################################
#                                                   #
# This file provides the numerical evidence for     #
# Conjecture 1 in the article. The details of the   #
# numerics are explained in Appendix B, section 2.  #
#                                                   #
# One can run this file by opening a terminal       #
# and typing                                        #
#                                                   #
#   python3 numerical_evidence_CKS2018_bounds.py    #
#                                                   #
# The results of the numerics will be stored in the #
# folder `output`.                                  #
#                                                   #
# This file is part of the supporting material      #
# belonging to:                                     #
# "Robust self-testing of two-qubit states"         #
# Tim Coopmans, Jędrzej Kaniewski and Christian     #
# Schaffner (2018)                                  #
# arXiv:                                            #
#                                                   #
#####################################################

import numpy as np
import CKS2018

def get_smallest_eigenvalue_of_hermitian(A):
    """
    Returns the smallest eigenvalue of 
    a hermitian matrix `A`.
    """
    (eigenvalues, __ ) = np.linalg.eig(A)
    smallest_eigenvalue = min(eigenvalues)

    # Since a hermitian matrix has real
    # eigenvalues, we may cast its 
    # smallest eigenvalue to a real number
    return np.real(smallest_eigenvalue)


def create_twodimensional_grid(x_start, x_stop, delta_x, y_start, y_stop, delta_y):
    number_of_x_grid_points = int(np.ceil( (x_stop - x_start) / delta_x))
    number_of_y_grid_points = int(np.ceil( (y_stop - y_start) / delta_y))
    x_points_list = np.linspace(start=x_start,
                                stop=x_stop,
                                num=number_of_x_grid_points,
                                endpoint=True)
    y_points_list = np.linspace(start=y_start,
                                stop=y_stop,
                                num=number_of_y_grid_points,
                                endpoint=True)
    return [(x, y) for x in x_points_list for y in y_points_list]

def get_minimum_of_smallest_eigenvalue_of_hermitian_over_grid(A, grid, **additional_parameters_for_A):
    """
    This method first iterates over points (a,b) 
    in a grid and at each point computes the smallest 
    eigenvalue of a hermitian matrix `A(a,b)`.
    It subsequently returns the minimum of
    these smallest eigenvalues.

    Parameters
    ----------
    A : method that takes parameters `a` and `b`, the grid
        coordinates, and possibly additional parameters
        as specified in `additional_parameters_for_A`.
    grid : list of tuples (`a`, `b`), where `a` and `b`
            are floats.
    """
    # set a start value
    a = grid[0][0]
    b = grid[0][1]
    current_minimum_eigenvalue = get_smallest_eigenvalue_of_hermitian(A(a=a, b=b, **additional_parameters_for_A))

    # compute the minimum of the smallest eigenvalue
    for (a, b) in grid:
        smallest = get_smallest_eigenvalue_of_hermitian(A(a=a, b=b, **additional_parameters_for_A))
        if current_minimum_eigenvalue > smallest:
            current_minimum_eigenvalue = smallest
    return current_minimum_eigenvalue

def write_data_to_csv_file(data_save_location, dictionary, description_preamble):
    with open(data_save_location, "w+") as datafile:
        datafile.write(description_preamble)
        for key, value in dictionary.items():
            datafile.write("{}, {}\n".format(key, value))
    datafile.close()

def get_description_string(epsilon, delta_alpha, worst_case_minimum_description):
    return \
"""#
# This file contains values for alpha in [0, 2), the parameter of
# the family of tilted CHSH operators W_{alpha}, and the minimum 
# eigenvalue of the operator 
#   T_{alpha} = K_{alpha} - s_{alpha} * W_{alpha} - mu_{alpha} * Identity,
# minimized over a discretization of [0, pi/2] x [0, pi/4].
# Here, K_{alpha}, s_{alpha} and mu_{alpha} are defined as in the 
# article.
#
# Used parameters:
# - epsilon : %s
# - delta_alpha : %s
#
# 
# The first column contains values for alpha, the second the minimum
# eigenvalue of T_{alpha} over the grid.
#%s#
"""%(epsilon, delta_alpha, worst_case_minimum_description)


if __name__ == "__main__":

    #################################
    # Parameters of the numerics    #
    #################################

    # Grid coarseness in both x and y
    # For the analysis of the article,
    # the parameter epsilon was set
    # to 0.01
    epsilon = 0.1 
    
    # For the analysis in the article,
    # the parameter delta_alpha
    # was set to 0.001
    delta_alpha = 0.1
    alpha_start = 0.
    alpha_stop = 2. - delta_alpha
    

    #################################
    # Performing the numerics       #
    #################################

    number_of_alphas = int(np.ceil( ( alpha_stop - alpha_start) / delta_alpha )) + 1
    alpha_list = np.linspace(start=alpha_start,
                             stop=alpha_stop,
                             num=number_of_alphas,
                             endpoint=True)
   

    # Create the grid: discretization of [0, pi/4] x [0, pi/2]
    grid = create_twodimensional_grid(x_start=0.,
                                      x_stop=np.pi/4.,
                                      delta_x = epsilon,
                                      y_start=0.,
                                      y_stop=np.pi/2.,
                                      delta_y = epsilon)


    # We store the minimum of the smallest eigenvalue for each alpha in a dictionary
    minevals = {}

    # Compute the minimum of the smallest eigenvalues for each alpha
    print("\nNow computing the minimum eigenvalues for:")
    for alpha in alpha_list:
        print("")
        print(r" +  \alpha = {} (will go up to but not including 2)".format(alpha))
        minevals[alpha] = get_minimum_of_smallest_eigenvalue_of_hermitian_over_grid(CKS2018.T,
                                                                                    grid=grid,
                                                                                    alpha=alpha)
        print("    Minimum of smallest eigenvalues: {}".format(minevals[alpha]))

    # Compute the "worst case", i.e. the minimum over all minima of smallest eigenvalues
    alpha_of_minimum = min(minevals, key=minevals.get)
    worst_case_minimum = minevals[alpha_of_minimum]
    worst_case_minimum_description=\
"""
# The minimum eigenvalue found over the entire grid is {}, which
# is achieved for alpha={}.
""".format(worst_case_minimum, alpha_of_minimum)
    print(worst_case_minimum_description)



    #################################
    # Storing the results           #
    #################################

    # Save the data
    description_preamble = get_description_string(epsilon=epsilon,
                                                  delta_alpha=delta_alpha,
                                                  worst_case_minimum_description=worst_case_minimum_description)
    data_save_location = "output/minimum_eigenvalues.csv"
    write_data_to_csv_file(data_save_location=data_save_location,
                           dictionary=minevals,
                           description_preamble=description_preamble)
    
    # User output
    print("The result of computing the minimum eigenvalues can be found in {}"\
           .format(data_save_location))

