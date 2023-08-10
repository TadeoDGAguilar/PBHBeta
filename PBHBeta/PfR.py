# PfR.py

from PBHBeta import functions
from PBHBeta import constants
import numpy as np

def get_P_k_RD(M, Beta, delta_c):
    """
    This functions obtain the Power Spectrum (PS) Value in standard Big Bang scenario.

    Parameters:

        M (numpy.ndarray): Corresponds to masses of PBHs, constrained and computed from any main functions.
        beta (numpy.ndarray): Corresponds to the abundance of PBHs obtained from any main functions.
        delta_c (float): Threshold amplitude for PBH formation.

    Returns:

        k (numpy.ndarray): Wave number values
        P_k (numpy.ndarray): Power Spectrum values

    Note:

        Before use this function, compute M and Beta for each standar scenario

    """
    k = functions.k_rad(np.array(M))
    sigma = np.array(functions.inverse_error(Beta, delta_c))
    P_k = constants.A * sigma ** 2
    return k, P_k