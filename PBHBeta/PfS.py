
from PBHBeta import functions
from PBHBeta import constants
from PBHBeta import BfS
import numpy as np

def k_rad(M):
    a_end_inf_rad = functions.a_endre(constants.rho_r0, constants.rho_end_inf)
    k_end = a_end_inf_rad*constants.H_end
    k_end_over_k_rad = (M/(7.1*10**-2*constants.gam_rad*(1.8*10**15/constants.H_end)))**(1/2)
    k = (k_end/k_end_over_k_rad)*constants.GeV*constants.metter_m1
    k = np.array(k)
    return k

def a_endinf(a_end_re,rho_end_re,rho_end_inf):
    return a_end_re*(rho_end_re/rho_end_inf)**(1./6)

def get_P_k_SD(M, N_stiff, omega, gam_sd):

    """
    This functions obtain the Power Spectrum (PS) constraint in a Stiff fluid dominated (SD) scenario

    Parameters:
        - M (numpy.ndarray): Corresponds to masses of PBHs, constrained and computed from any main functions.
        - N_stiff (float): Is the total number of e-folds that the stiff era lasted.
        - omega (float): This value is to assign the equation of state
        - gam_sd: The particular value of \gamma^{SD} is not well known and we thus adopt \gamma^{SD} = 1

    Returns:

        - k (numpy.ndarray): Wave number values in SD.
        - P_k (numpy.ndarray): Power Spectrum values in SD.
        - betas_stiff (numpy.ndarray): Abundances of PBHs in SD.

    """

    C = (8/2)**(2/4)
    A = (5+3/3)**2/(4*(1+1/3)**2)
    betas_stiff = BfS.get_betas_stiff_tot(M, N_stiff, omega, gam_sd)
    betas_full = functions.get_Betas_full(M)
    k_end_over_k_stiff = (M/(7.1*10**-2*gam_sd*(1.8*10**15/constants.H_end)))**(2/3)
    sigma_tot = np.array(functions.inverse_error(betas_full, 0.41))
    a_end_reh = functions.a_endre(constants.rho_r0, constants.rho_end_inf * np.exp(-6*N_stiff))
    a_end_inf = a_endinf(a_end_reh, constants.rho_end_inf * np.exp(-6*N_stiff), constants.rho_end_inf)
    k_phys_rad = np.array(k_rad(M))
    sigma_st = M*0
    k_st = M*0
    for i in range(len(M)):
        if betas_stiff[i] != betas_full[i]:
            k_st[i] = a_end_inf*constants.H_end*constants.GeV*constants.metter_m1/k_end_over_k_stiff[i]
            sigma_st[i] = np.array(functions.inverse_error([betas_stiff[i]], 0.52))
        else:
            k_st[i] = k_phys_rad[i]
            sigma_st[i] = (A/C)**(1/2)*sigma_tot[i]

    P_k_st = constants.C*sigma_st**2

    return k_st, P_k_st, betas_stiff