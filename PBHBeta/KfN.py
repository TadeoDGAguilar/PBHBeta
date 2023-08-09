#Module KfN.py
#This module permit to the user compute to get modes K from Numer of efolds

from PBHBeta import functions
from PBHBeta import constants
from PBHBeta import BfN
import numpy as np
from scipy.optimize import fsolve


def k_MD(M):
  a_end_inf_rad = functions.a_endre(constants.rho_r0, constants.rho_end_inf)
  k_end = a_end_inf_rad * constants.H_end
  k_end_over_k_rad = (M / (7.1 * 10 ** -2 * constants.gam_rad * (1.8 * 10 ** 15 / constants.H_end))) ** (1 / 2)
  k = (k_end / k_end_over_k_rad) * constants.GeV * constants.metter_m1
  k = np.array(k)
  return k

def a_endinf(a_end_re, rho_end_re, rho_end_inf):
  return a_end_re * (rho_end_re / rho_end_inf) ** (1. / 3)

def get_P_k_MD(M, N_re, omega, gam_reh):


  def equation(p):
    x = p[0]
    return[1.9*10**-7*x**2*np.exp(-0.1*x**(-2/3))-beta_t]

  rho_end_inf = constants.rho_end_inf
  k_md = np.array(k_MD(M))
  betas_reh = BfN.get_betas_reh_tot(N_re, omega, gam_reh)
  k_end_over_k_reh = (M/(7.1*10**-2*gam_reh*(1.8*10**15/constants.H_end)))**(1/3)
  betas_full = functions.get_Betas_full(M)
  sigma_tot = np.array(functions.inverse_error(betas_full, 0.41))
  k = M*0
  sigma = M*0
  a_end_reh = functions.a_endre(constants.rho_r0,constants.rho_end_inf*np.exp(-3*N_re))
  a_end_inf = a_endinf(a_end_reh,constants.rho_end_inf*np.exp(-3*N_re),rho_end_inf)
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
  return sigma, k, betas_reh