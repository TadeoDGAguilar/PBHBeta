# File: BfS.py
## Module with BfS to compute
### Version BETA

from scipy.optimize import fsolve
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
from numpy import diff
import sys
from scipy import special
from PBHBeta import constants
from PBHBeta import constraints
from PBHBeta.functions import get_Betas_full
from PBHBeta.functions import get_Omegas_full
import numpy as np

# Set the values of the variables
rho_end_inf =  constants.rho_end_inf #define rho_end_inf
H_end = constants.H_end #define H_end
gam_rad = constants.gam_rad
rho_end = constants.rho_end # define rho_end
M_pl_g = constants.M_pl_g # define M_pl_g
ln_den_end = np.log(constants.rho_end)
t_pl_s = constants.t_pl_s # s
s_to_evm1 = constants.s_to_evm1
t_pl = constants.t_pl
M_pl = constants.M_pl
ev1 = constants.ev1
ev2 = constants.ev2
rho_end = constants.rho_end
log10_M_tot = np.linspace(0,20,1000)
ln_den_end = np.log(rho_end)


def get_betas_stiff_tot(M, N_stiff, omega, gam_stiff):
    """
    This function contains compute the abundance of PBHs in a Stiff-fluid Domination (SD) era. This includes the system of equations that describe the evolution of the population of PBHs after their formation.
    Eqs. (35) to (40) from the article.

    Parameters:
        - M (array): Range of Mass from PBHs
        - N_stiff (float): Is the total number of e-folds that the stiff era lasted.
        - omega (float): This value is to assign the equation of state
        - gam_sd: It is a constant that encrypts the efficiency of the collapse and for a Stiff Fluid Domination era. The particular value of \gamma^{SD} is not well known and we thus adopt \gamma^{SD} = 1

    Returns:
        - betas_stiff (list): Contain the total abundances of PBHs in SD.

    """
    M_tot = np.array(M)
    Omegas_full = get_Omegas_full(M)#np.array(constraints.Omegas_full)
    betas_full = get_Betas_full(M)#np.array(constraints.betas_full)
    
    rho_end_stiff = rho_end_inf*np.exp(-6*N_stiff)
    if rho_end_stiff <= rho_end:
        raise ValueError("The end of the stiff era happens after BBN")
    ln_den_end_stiff = np.log(rho_end_stiff)
    k_end_over_k_stiff = (M_tot/(7.1*10**-2*gam_stiff*(1.8*10**15/H_end)))**(2/3)
    rho_form_stiff = rho_end_inf/(k_end_over_k_stiff)**3
    ln_rho_form_stiff = np.log(rho_form_stiff)
    
    betas_stiff = []
    
    for i in range(len(M_tot)):
        end_value = Omegas_full[i]
        M = M_tot[i]
        ln_den_f = np.log(rho_form_stiff[i])
        
        def diff_ext_rel(ln_rho,initial,M,beta0,omega):
            dy = np.zeros(initial.shape)
            b = initial[0]
            Om_ext = initial[1]
            Om_0 = beta0*b*(M_pl_g/M)
            dy[0] = -(Om_0+(1-3*omega)*Om_ext-1.)*b/(Om_0+(1-3*omega)*Om_ext-4.)
            dy[1] = -(Om_0+(1-3*omega)*Om_ext+3*omega-1)*Om_ext/(Om_0+(1-3*omega)*Om_ext-4)
            return dy
        
        def diff_ext(ln_rho,initial,M,beta0,omega):
            dy = np.zeros(initial.shape)
            b = initial[0]
            time = initial[1]
            Om_ext = initial[2]
            Delta_t = t_pl*(M/M_pl_g)**3
            Om_0 =  beta0*b*(1.-time/Delta_t)**(1./3)
            dy[0] = -(Om_0+(1-3*omega)*Om_ext-1.)*b/(Om_0+(1-3*omega)*Om_ext-4.)
            dy[1] = 3**(1./2)*M_pl/((Om_0+(1-3*omega)*Om_ext-4.)*np.exp(ln_rho)**(1./2))
            dy[2] = -(Om_0+(1-3*omega)*Om_ext+3*omega-1)*Om_ext/(Om_0+(1-3*omega)*Om_ext-4)
            return dy
        
        def end_evol_ext(ln_rho,initial,M,beta0,omega):
            Delta_t = t_pl*(M/M_pl_g)**3
            Mass_end = M**(1.-diff_ext(ln_rho,initial,M,beta0,omega)[1]/Delta_t)**(1./3)
            return Mass_end - M_pl_g
        
        def objective_stiff(beta0):
            N_f = (1/6)*(ln_den_f-ln_den_end_stiff)
            initial_stiff = np.array([1.,0.,(1.-beta0[0]*(1+np.exp(N_f)))/(1+np.exp(-2*N_f))])
            if initial_stiff[1]<1/2:
                initial_stiff[1] = 1/2
            sol = solve_ivp(diff_ext,(ln_den_f,ln_den_end),initial_stiff,t_eval=ln_den,args=(M,beta0[0],omega),method = "DOP853")
            Delta_t = t_pl*(M/M_pl_g)**3
            y = beta0[0]*sol.y[0]*(1.-sol.y[1]/Delta_t)**(1./3)
            return y[-1]-end_value
        
        def objective_stiff_rel(beta0):
            N_f = (1/6)*(ln_den_f-ln_den_end_stiff)
            initial_stiff = np.array([1.,(1.-(M_pl_g/M)*beta0[0]*(1+np.exp(N_f)))/(1+np.exp(-2*N_f))])
            sol = solve_ivp(diff_ext_rel,(ln_den_f,ln_den_end),initial_stiff,t_eval=ln_den,args=(M,beta0[0],omega),method = "DOP853")#
            y = beta0[0]*sol.y[0]*(M_pl_g/M)
            return y[-1]-end_value
        
        if ln_den_f <= ln_den_end_stiff:
            betas_stiff.append(betas_full[i])
            continue
        ln_den = np.linspace(ln_den_f, ln_den_end, 10000)
        sol_try = solve_ivp(diff_ext, (ln_den_f, ln_den_end), np.array([1., 0., (1-1e-10)/(1+np.exp(-2*N_stiff)/0.1)]), events=end_evol_ext, t_eval=ln_den, args=(M, 1e-10, omega), method="DOP853")
        if sol_try.t[-1] > ln_den_end:
            beta_try = 2*Omegas_full[i]*(M/M_pl_g)/(1+np.sqrt(rho_form_stiff[i]/(rho_end_stiff))*(rho_end_stiff/rho_end)**(1/4))
            beta0, = fsolve(objective_stiff_rel, beta_try)
            if (beta0/beta_try > 10000):
                betas_stiff.append(beta_try)
            else:
                betas_stiff.append(beta0)
        else:
            beta_try = 1e-4*Omegas_full[i]/(1+np.sqrt(rho_form_stiff[i]/(rho_end_stiff))*(rho_end_stiff/rho_end)**(1/4))
            beta0, = fsolve(objective_stiff, beta_try)
            betas_stiff.append(beta0)
    return betas_stiff




