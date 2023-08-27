# PfR.py

from PBHBeta import functions
from PBHBeta import constants
import numpy as np

def k_RD(M):
    a_end_inf_rad = (constants.rho_r0 / constants.rho_end_inf) ** (1. / 4)
    k_end = a_end_inf_rad * constants.H_end
    k_end_over_k_rad = (M/(7.1*10**-2*constants.gam_rad*(1.8*10**15/constants.H_end)))**(1/2)
    k = (k_end/k_end_over_k_rad)*constants.GeV*constants.metter_m1
    k = np.array(k)
    return k

def get_P_k_RD(M, Beta, delta_c):
    """
    This functions obtain the Power Spectrum (PS) Value in standard Big Bang scenario.

    Parameters:

        - M (numpy.ndarray): Corresponds to masses of PBHs, constrained and computed from any main functions.
        - beta (numpy.ndarray): Corresponds to the abundance of PBHs obtained from any main functions.
        - delta_c (float): Threshold amplitude for PBH formation.

    Returns:

        - k (numpy.ndarray): Wave number values
        - P_k (numpy.ndarray): Power Spectrum values

    Note:

        Before use this function, compute M and Beta for each standar scenario

    """
    k = k_RD(np.array(M))
    sigma = np.array(functions.inverse_error(Beta, delta_c))
    P_k = constants.A * sigma ** 2
    return k, P_k