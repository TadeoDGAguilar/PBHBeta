#Module PfM.py
#This module permit to the user compute to get modes K from Numer of efolds

from PBHBeta import functions
from PBHBeta import constants
from PBHBeta import BfM
import numpy as np
from scipy.optimize import fsolve


def k_MD(M):
  a_end_inf_rad = functions.a_endre(constants.rho_r0, constants.rho_end_inf)
  k_end = a_end_inf_rad * constants.H_end
  k_end_over_k_rad = (M / (7.1 * 10 ** -2 * constants.gam_rad * (1.8 * 10 ** 15 / constants.H_end))) ** (1 / 2)
  k = (k_end / k_end_over_k_rad) * constants.GeV * constants.metter_m1
  k = np.array(k)
  return k

def a_endre(rho_r0, rho_end_re):
  return (rho_r0 / rho_end_re) ** (1. / 4)

def a_endinf(a_end_re, rho_end_re, rho_end_inf):
  return a_end_re * (rho_end_re / rho_end_inf) ** (1. / 3)

def get_P_k_MD(M, N_md, omega, gam_md):

  """
  This functions obtain the Power Spectrum (PS) constraint in a Early matter dominated (MD) scenario

  Parameters:

      - M (numpy.ndarray): Corresponds to masses of PBHs, constrained and computed from any main functions.
      - N_md (float): This is the total number of e−folds of MD.
      - omega (float): This value is to assign the equation of state
      - gam_md: . The particular value of \gamma^{MD} is not well known and we thus adopt \gamma^{MD} = 1

  Returns:

      - k (numpy.ndarray): Wave number values in MD.
      - P_k (numpy.ndarray): Power Spectrum values in MD.
      - betas_reh (numpy.ndarray): Abundances of PBHs in MD.

  """

  def equation(p):
    x = p[0]
    return[1.9*10**-7*x**2*np.exp(-0.1*x**(-2/3))-beta_t]

  rho_end_inf = constants.rho_end_inf
  k_md = np.array(k_MD(M))
  betas_reh = BfM.get_betas_reh_tot(M, N_md, omega, gam_md)
  k_end_over_k_reh = (M/(7.1*10**-2*gam_md*(1.8*10**15/constants.H_end)))**(1/3)
  betas_full = functions.get_Betas_full(M)
  sigma_tot = np.array(functions.inverse_error(betas_full, 0.41))
  k = M*0
  sigma = M*0
  a_end_reh = functions.a_endre(constants.rho_r0,constants.rho_end_inf*np.exp(-3*N_md))
  a_end_inf = a_endinf(a_end_reh,constants.rho_end_inf*np.exp(-3*N_md),rho_end_inf)
  for i in range(len(M)):
    if betas_reh[i] != betas_full[i]:
      k[i] = a_end_inf*constants.H_end*constants.GeV*constants.metter_m1/k_end_over_k_reh[i]
      sigma_try = (betas_reh[i]/0.05556)**(1/5)
      if sigma_try >= 0.005:
        sigma[i] = sigma_try
      else:
        beta_t = betas_reh[i]
        sigma[i] = fsolve(equation,[0.005])
    else:
        k[i] = k_md[i]
        sigma[i] = sigma_tot[i]

  P_k_md = constants.B*sigma**2
  return k, P_k_md, betas_reh