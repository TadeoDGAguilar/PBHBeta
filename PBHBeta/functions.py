# File: functions.py
## Module with functions
### Version BETA

from scipy.optimize import fsolve
from scipy.optimize import least_squares
from scipy.integrate import solve_ivp
from numpy import diff
import sys
from scipy import special
from PBHBeta import constants
from PBHBeta import constraints
import numpy as np

def put_M_array(Mass_min, Mass_max):
    """
    Generate an array of primordial black hole (PBH) masses in grams based on specified limits.

    Parameters:
        Mass_min (float): The minimum PBH mass value in grams.
        Mass_max (float): The maximum PBH mass value in grams.

    Returns:
        np.ndarray: An array of calculated PBH masses.
    """
    i = 0
    M = 0
    delta_M = 0.0123
    M_tot_try = []
    num_values = 20

    mass_array = np.geomspace(Mass_min, 10**(i*delta_M) , num_values)

    while M < constraints.data_mass[0]:
        M = 10**(i*delta_M)
        M_tot_try.append(M)
        i = i+1

    M_tot_try = np.concatenate((mass_array, M_tot_try, constraints.data_mass))
    M_tot_try = np.unique(M_tot_try)
    M = M_tot_try[-1]

    A = M
    j = 0

    while M < Mass_max:
        j = j+1
        M = A*10**(j*delta_M)
        M_tot_try = np.append(M_tot_try,[M])

    constraints.M_tot = np.array(M_tot_try)

    return constraints.M_tot


def diff_rad_rel(ln_rho,initial,M,beta0):

    """In the scenario where PBHs evaporate before reaching the energy scale of interest (as is the case, for example,
    before reaching the energy scale of BBN), we calculate the PBH abundance by assuming the existence of remnants with
    a mass equal to the Planck mass. Instead of simultaneously solving Eqs.(10) and (11) with the constraint Eq. (8), we
    focus on solving Eq.(10) with the constraint $\Omega_{PBH} = (m_{Pl}/M_{PBH})\beta(M_{PBH})$."""

    # Extract initial scale factor b and calculate Om_0
    b = initial[0]
    Om_0 = beta0 * b * (constants.M_pl_g / M)

    # Calculate the derivative of the scale factor b
    dy = -(Om_0 - 1.) * b / (Om_0 - 4.)

    return dy


def diff_rad(ln_rho,initial,M,beta0):
    """This function corresponds to Eqs.(10) and (11) with the constraint Eq.(8) in our reference paper. It is employed
    to calculate the abundance of PBHs in a radiation-dominated universe as a function of total energy density."""

    # Initialize dy array
    dy = np.zeros(initial.shape)

    # Extract initial values of scale factor b and time
    b = initial[0]
    time = initial[1]

    # Calculate Delta_t and Om_0
    Delta_t = constants.t_pl * (M / constants.M_pl_g) ** 3
    Om_0 = beta0 * b * (1. - time / Delta_t) ** (1. / 3)

    # Calculate the derivative of the scale factor b and the time derivative of the density of radiation
    dy[0] = -(Om_0 - 1.) * b / (Om_0 - 4.)
    dy[1] = 3 ** (1. / 2) * constants.M_pl / ((Om_0 - 4.) * np.exp(ln_rho) ** (1. / 2))

    return dy


def end_evol(ln_rho,initial,M,beta0):
    """This function is used to determine whether a PBH reaches the Planck mass (thus becoming a Planck relic) or not.
    By solving the system of equations (10) and (11) with the constraint (8) from our reference article, this function
    is used as a stopping condition for the evolution of the system. In the event that the evolution is halted before
    reaching the desired energy scale (such as the scale of BBN), the evolution of PBHs is carried out considering them
    as Planck mass relics."""
    # Calculate Delta_t and Mass_end
    Delta_t = constants.t_pl * (M / constants.M_pl_g) ** 3
    Mass_end = M * (1. - diff_rad(ln_rho,initial,M,beta0)[1] / Delta_t) ** (1. / 3)

    # Return the difference between the final mass of a system and the Planck mass
    return Mass_end - constants.M_pl_g



def k_end_over_k(Mpbh, omega):
    """
    Calculates the ratio of k_end/k for a given PBH mass and radiation energy density parameter.

    Parameters:
        - Mpbh (float): The mass of the PBH, in grams.
        - omega (float): The energy density parameter for radiation.

    Returns:
        - ratio (float): The ratio of k_end/k for the given PBH mass and radiation energy density parameter.
    """
    if omega==1/3:
        res = (Mpbh/(7.1*10**-2*constants.gam_rad*(1.8*10**15/constants.H_end)))**(1/2)
    else:
        z = (1+3*omega)/(3*(1 + omega))
        ratio = (Mpbh*constants.H_end/(3*constants.gam_rad*(constants.M_pl**2.)))**z
        res = np.array(ratio)
    return res


def rho_f(Mpbh, omega):
    """
    Calculates the final density of black holes after evaporation.

    Parameters:
        - Mpbh (float): The initial mass of a black hole, in grams.
        - omega (float): The ratio of the energy density of dark matter to the critical density of the universe.

    Returns:
        - rho (float): The final density of black holes, in grams per cubic centimeter.
    """
    if omega==1/3:
        k_end_over_k_rad = (Mpbh/(7.1*10**-2*constants.gam_rad*(1.8*10**15/constants.H_end)))**(1/2)
        rho_f = constants.rho_end_inf/(k_end_over_k_rad)**4
    else:
        z = (1+3*omega)/(3*(1 + omega))
        ratio = (Mpbh*constants.H_end/(3*constants.gam_rad*(constants.M_pl**2.)))**z
        res = np.array(ratio)
        i = (6*(1 + omega))/(1+(3*omega))
        rho_f = constants.rho_end_inf/(res)**i
    return rho_f



ln_den_end = np.log(constants.rho_end)


def Betas_DM(M_tot, omega):
    """
    Calculates the abundance of PBHs for dark matter constraints. See equations (13) and (19).

    Parameters:
        - M_tot (array-like): Array of masses in grams.
        - omega (float): This value is to assign the equation of state


    Returns:
        A tuple containing four numpy arrays:
            - M_n (numpy.ndarray): Represents the masses of PBHs can be considered candidates for Dark Matter (DM).
            - beta (numpy.ndarray): Corresponds to the abundance obtained from M_n.
            - M_relic (numpy.ndarray): Masses of the relic dark matter components.
            - beta_relic_prim (numpy.ndarray): Corresponds to abundance obtained from M_relic.
            - Omegas_tot (numpy.ndarray): Correspond to evolution of abundance of PBHs (with DM constraint) after their formation.
    """

    M_n = []
    betas_prim = []
    M_relic = []
    betas_relic_prim = []
    betas_tot = []
    Omegas_tot = []
    Omegas = []
    Omegas_relic_pbbn = []
    Omegas_relic = []
    M_dm = []
    M_dm_rel_pbbn = []
    M_dm_rel = []

    M_pl = 1.22089 * 10 ** 19
    M_pl_g = 2.17645e-5
    t_pl_s = 5.39 * 10 ** -44
    s_to_evm1 = (1. / 6.5823) * 10 ** 25
    t_pl = t_pl_s * s_to_evm1

    gam_rad = (1. / 3) ** (3. / 2)
    H_end = 4.44 * 10 ** 13.
    rho_end_inf = 3. * M_pl ** 2. * H_end ** 2

    rho_form_rad = rho_f(M_tot, omega)

    # k_end_over_k_rad = (M_tot * H_end / (gam_rad * 3 * M_pl))**(1 + 3 * w) / (3 * (1 + w))
    # rho_form_rad = rho_end_inf / (k_end_over_k_rad)**((6 * (1 + w)) / (1 + (3 * w)))
    rho_end = (1e-2) ** 4
    ln_den_end = np.log(rho_end)
    rho_end = (1e-2) ** 4
    ln_den_end = np.log(rho_end)

    for i in range(len(M_tot)):
        if M_tot[i] > 4.1 * 10 ** 14:
            M_n.append(M_tot[i])
            beta = 1.86 * 10 ** -18 * (M_tot[i] / (10 ** 15)) ** (1 / 2)
            betas_prim.append(beta)
        elif M_tot[i] < 10 ** 11 * constants.M_pl_g:
            M_relic.append(M_tot[i])
            beta = 2 * 10 ** -28 * (M_tot[i] / constants.M_pl_g) ** (3 / 2)
            betas_relic_prim.append(beta)
        else:
            beta = constants.ev1
        betas_tot.append(beta / gam_rad ** (1 / 2))
    betas_prim = np.array(betas_prim)
    betas_relic_prim = np.array(betas_relic_prim)
    betas_tot = np.array(betas_tot)
    constraints.betas_DM_tot = betas_tot

    M_n = np.array(M_n)
    M_relic = np.array(M_relic)
    betas = betas_prim / constants.gam_rad ** (1 / 2)
    betas_relic = betas_relic_prim / constants.gam_rad ** (1 / 2)
    for i in range(len(constraints.betas_DM_tot)):
        ln_den_f = np.log(rho_form_rad[i])
        if ln_den_f <= ln_den_end:
            continue
        ln_den = np.linspace(ln_den_f, ln_den_end, 10000)
        sol_try = solve_ivp(diff_rad, (ln_den_f, ln_den_end), np.array([1., 0.]), events=end_evol, t_eval=ln_den,
                            args=(M_tot[i], betas_tot[i]), method="DOP853")
        if sol_try.t[-1] > ln_den_end:
            sol_try = solve_ivp(diff_rad_rel, (ln_den_f, ln_den_end), np.array([1.]), t_eval=ln_den,
                                args=(M_tot[i], betas_tot[i]), method="DOP853")
            y = betas_tot[i] * sol_try.y[0][-1] * (M_pl_g / M_tot[i])
            if M_tot[i] < 10 ** 11 * M_pl_g:
                Omegas_relic_pbbn.append(y)
                M_dm_rel_pbbn.append(M_tot[i])
        else:
            Delta_t = t_pl * (M_tot[i] / M_pl_g) ** 3
            y = betas_tot[i] * sol_try.y[0][-1] * (1. - sol_try.y[1][-1] / Delta_t) ** (1. / 3)
            if M_tot[i] > 4.1 * 10 ** 14:
                Omegas.append(y)
                M_dm.append(M_tot[i])
            elif M_tot[i] < 10 ** 11 * M_pl_g:
                Omegas_relic.append(y)
                M_dm_rel.append(M_tot[i])
        Omegas_tot.append(y)

    Omegas_tot = np.array(Omegas_tot)
    constraints.Omega_DM_tot = Omegas_tot
    Omegas = np.array(Omegas)
    Omegas_relic_pbbn = np.array(Omegas_relic_pbbn)
    Omegas_relic = np.array(Omegas_relic)
    M_dm = np.array(M_dm)
    M_dm_rel_pbbn = np.array(M_dm_rel_pbbn)
    M_dm_rel = np.array(M_dm_rel)

    return M_n, betas, M_relic, betas_relic, Omegas_tot



    
def Betas_BBN(M_tot, omega):
    """
    Calculates the abundance of PBHs for Big Bang Nucleosynthesis constraints.

    Parameters:
        - M_tot (array-like): Array of masses in grams.
        - omega (float): This value is to assign the equation of state.

    Returns:
        A tuple containing three numpy arrays:
            - M_bbn (numpy.ndarray): Represents the masses of PBHs leading to an injection of particles that can impact the predictions of BBN.
            - betas_bbn (numpy.ndarray): Corresponds to the abundance obtained from M_bbn.
            - Omegas_bbn_tot (numpy.ndarray): Corresponds to evolution of abundance of PBHs with BBN constraint.
    """

    betas_bbn = []
    M_bbn = []
    M_bbn_bbn = []
    Omegas_bbn = []
    Omegas_bbn_tot = []
    Omegas_bbn_pbbn = []
    M_bbn_pbbn = []

    rho_form_rad = rho_f(M_tot, omega)
    j = 0
    
    for i in range(len(M_tot)):
        if M_tot[i]>= constraints.data_mass[0] and M_tot[i]<constraints.data_mass[76]:
            M_bbn.append(M_tot[i])
            beta = constraints.data_abundances[j]/constants.gam_rad**(1/2)
            betas_bbn.append(beta)
        
            ln_den_f = np.log(rho_form_rad[i])
            if ln_den_f <= ln_den_end:
                continue
            ln_den = np.linspace(ln_den_f,ln_den_end,10000)
            sol_try = solve_ivp(diff_rad,(ln_den_f,ln_den_end),np.array([1.,0.]),events=end_evol,t_eval=ln_den,args=(M_tot[i],beta),method = "DOP853")
            if sol_try.t[-1] > ln_den_end:
                sol_try = solve_ivp(diff_rad_rel,(ln_den_f,ln_den_end),np.array([1.]),t_eval=ln_den,args=(M_tot[i],beta),method = "DOP853")
                y = beta*sol_try.y[0][-1]*(constants.M_pl_g/M_tot[i])
                Omegas_bbn_pbbn.append(y)
                M_bbn_pbbn.append(M_tot[i])
            else:
                Delta_t = constants.t_pl*(M_tot[i]/constants.M_pl_g)**3
                y = beta*sol_try.y[0][-1]*(1.-sol_try.y[1][-1]/Delta_t)**(1./3)
                Omegas_bbn.append(y)
                M_bbn_bbn.append(M_tot[i])
            j = j+1
        elif M_tot[i]>= constraints.data_mass[76] and M_tot[i]<2.5*10**13:
            M_bbn.append(M_tot[i])
            beta = constraints.data_abundances[76]/constants.gam_rad**(1/2)
            betas_bbn.append(beta)
        
            ln_den_f = np.log(rho_form_rad[i])
            if ln_den_f <= ln_den_end:
                continue
            ln_den = np.linspace(ln_den_f,ln_den_end,10000)
            sol_try = solve_ivp(diff_rad,(ln_den_f,ln_den_end),np.array([1.,0.]),events=end_evol,t_eval=ln_den,args=(M_tot[i],beta),method = "DOP853")
            Delta_t = constants.t_pl*(M_tot[i]/constants.M_pl_g)**3
            y = beta*sol_try.y[0][-1]*(1.-sol_try.y[1][-1]/Delta_t)**(1./3)
            Omegas_bbn.append(y)
            M_bbn_bbn.append(M_tot[i])
            j = j+1
        else:
            beta = constants.ev1
            y = constants.ev2
        constraints.betas_BBN_tot.append(beta)
        Omegas_bbn_tot.append(y)

    betas_bbn = np.array(betas_bbn)
    Omegas_bbn_tot = np.array(Omegas_bbn_tot)
    M_bbn = np.array(M_bbn)
    M_bbn_bbn = np.array(M_bbn_bbn)
    M_bbn_pbbn = np.array(M_bbn_pbbn)
    Omegas_bbn = np.array(Omegas_bbn)
    Omegas_bbn_pbbn = np.array(Omegas_bbn_pbbn)
    constraints.Omega_BBN_tot = Omegas_bbn_tot
    
    return M_bbn, betas_bbn, Omegas_bbn_tot



def Betas_SD(M_tot, omega):

    """
    Calculates the abundance of PBHs for Spectral Distortion (SD) constraints.

    Parameters:
        - M_tot (array-like): Array of masses in grams.
        - omega (float): This value is to assign the equation of state.

    Returns:
        A tuple containing three numpy arrays:
            - M_sd (numpy.ndarray): Represents the masses of PBHs consider the constraints from SD
            - betas_sd (numpy.ndarray): Corresponds to the abundance obtained from M_sd.
            - Omegas_sd (numpy.ndarray): Corresponds to evolution of abundance of PBHs with SD constraint.
    """

    
    betas_sd = []
    M_sd = []
    M_sd_bbn = []
    Omegas_sd = []
    
    rho_form_rad = rho_f(M_tot, omega)

    for i in range(len(M_tot)):
        if M_tot[i]> 10**11 and M_tot[i]<10**13:
            M_sd.append(M_tot[i])
            beta = 10**(-21)/constants.gam_rad**(1/2)
            betas_sd.append(beta)
        
            ln_den_f = np.log(rho_form_rad[i])
            if ln_den_f <= ln_den_end:
                continue
            ln_den = np.linspace(ln_den_f,ln_den_end,10000)
            sol_try = solve_ivp(diff_rad,(ln_den_f,ln_den_end),np.array([1.,0.]),events=end_evol,t_eval=ln_den,args=(M_tot[i],beta),method = "DOP853")
            Delta_t = constants.t_pl*(M_tot[i]/constants.M_pl_g)**3
            y = beta*sol_try.y[0][-1]*(1.-sol_try.y[1][-1]/Delta_t)**(1./3)
            Omegas_sd.append(y)
            M_sd_bbn.append(M_tot[i])
        else:
            beta = constants.ev1
            y = constants.ev2
        constraints.betas_SD_tot.append(beta)
        constraints.Omega_SD_tot.append(y)
    
    betas_sd = np.array(betas_sd)
    M_sd = np.array(M_sd)
    M_sd_bbn = np.array(M_sd_bbn)
    Omegas_sd = np.array(Omegas_sd)
    
    return M_sd, betas_sd, Omegas_sd



def Betas_CMB_AN(M_tot, omega):

    """
    Calculates the abundance of PBHs for CMB Anisotropies (CMB_AN) constraints.

    Parameters:
            - M_tot (array-like): Array of masses in grams.
            - omega (float): This value is to assign the equation of state.

    Returns:
        A tuple containing three numpy arrays:
            - M_an (numpy.ndarray): Represents the masses of PBHs consider the constraints from CMB_AN.
            - betas_an (numpy.ndarray): Corresponds to the abundance obtained from M_an.
            - Omegas_an (numpy.ndarray): Corresponds to evolution of abundance of PBHs with CMB_AN constraint.
    """


    betas_an = []
    M_an = []
    M_an_bbn = []
    Omegas_an = []
    
    rho_form_rad = rho_f(M_tot, omega)

    for i in range(len(M_tot)):
        if M_tot[i]> 2.5*10**13 and M_tot[i]<2.4*10**14:
            M_an.append(M_tot[i])
            beta = 3*10**(-30)*(M_tot[i]/10**13)**3.1/constants.gam_rad**(1/2)
            betas_an.append(beta)
        
            ln_den_f = np.log(rho_form_rad[i])
            if ln_den_f <= ln_den_end:
                continue
            ln_den = np.linspace(ln_den_f,ln_den_end,10000)
            sol_try = solve_ivp(diff_rad,(ln_den_f,ln_den_end),np.array([1.,0.]),events=end_evol,t_eval=ln_den,args=(M_tot[i],beta),method = "DOP853")
            Delta_t = constants.t_pl*(M_tot[i]/constants.M_pl_g)**3
            y = beta*sol_try.y[0][-1]*(1.-sol_try.y[1][-1]/Delta_t)**(1./3)
            Omegas_an.append(y)
            M_an_bbn.append(M_tot[i])
        else:
            beta = constants.ev1
            y = constants.ev2
        constraints.betas_CMB_AN_tot.append(beta)
        constraints.Omega_CMB_AN_tot.append(y)
    
    betas_an = np.array(betas_an)
    M_an = np.array(M_an)
    M_an_bbn = np.array(M_an_bbn)
    Omegas_an = np.array(Omegas_an)
    
    return M_an, betas_an, Omegas_an


def Betas_GRB(M_tot, omega):

    """
    Calculates the abundance of PBHs could contribute to the diffuse x-ray and γ-ray background (GRB).

        Parameters:
                - M_tot (array-like): Array of masses in grams.
                - omega (float): This value is to assign the equation of state.

        Returns:
            A tuple containing six numpy arrays:
                - M_grb1 (numpy.ndarray): Represents the masses of PBHs in the range [3x10^(13), 5.1x10^(14)].
                - M_grb2 (numpy.ndarray): Represents the masses of PBHs in the range [5.1x10^(14), 7x10^(16)].
                - betas_grb1 (numpy.ndarray): Corresponds to the abundance obtained from M_grb1.
                - betas_grb2 (numpy.ndarray): Corresponds to the abundance obtained from M_grb2.
                - Omegas_grb1 (numpy.ndarray): Corresponds to evolution of abundance of betas_grb1.
                - Omegas_grb2 (numpy.ndarray): Corresponds to evolution of abundance of betas_grb2.
    """

    betas_grb1 = []
    M_grb1 = []
    betas_grb2 = []
    M_grb2 = []

    M_grb1_bbn = []
    Omegas_grb1 = []
    M_grb2_bbn = []
    Omegas_grb2 = []
    
    rho_form_rad = rho_f(M_tot, omega)

    for i in range(len(M_tot)):
        if M_tot[i]> 3*10**13 and M_tot[i]<4.1*10**14:
            M_grb1.append(M_tot[i])
            beta = 5*10**(-28)*(M_tot[i]/(4.1*10**14))**-3.3/constants.gam_rad**(1/2)
            betas_grb1.append(beta)
        
            ln_den_f = np.log(rho_form_rad[i])
            if ln_den_f <= ln_den_end:
                continue
            ln_den = np.linspace(ln_den_f,ln_den_end,10000)
            sol_try = solve_ivp(diff_rad,(ln_den_f,ln_den_end),np.array([1.,0.]),events=end_evol,t_eval=ln_den,args=(M_tot[i],beta),method = "DOP853")
            Delta_t = constants.t_pl*(M_tot[i]/constants.M_pl_g)**3
            y = beta*sol_try.y[0][-1]*(1.-sol_try.y[1][-1]/Delta_t)**(1./3)
            Omegas_grb1.append(y)
            M_grb1_bbn.append(M_tot[i])
        
        elif M_tot[i]> 4.1*10**14 and M_tot[i]<7*10**16:
            M_grb2.append(M_tot[i])
            beta = 5*10**(-26)*(M_tot[i]/(4.1*10**14))**3.9/constants.gam_rad**(1/2)
            betas_grb2.append(beta)
        
            ln_den_f = np.log(rho_form_rad[i])
            if ln_den_f <= ln_den_end:
                continue
            ln_den = np.linspace(ln_den_f,ln_den_end,10000)
            sol_try = solve_ivp(diff_rad,(ln_den_f,ln_den_end),np.array([1.,0.]),events=end_evol,t_eval=ln_den,args=(M_tot[i],beta),method = "DOP853")
            Delta_t = constants.t_pl*(M_tot[i]/constants.M_pl_g)**3
            y = beta*sol_try.y[0][-1]*(1.-sol_try.y[1][-1]/Delta_t)**(1./3)
            Omegas_grb2.append(y)
            M_grb2_bbn.append(M_tot[i])
        else:
            beta = constants.ev1
            y = constants.ev2
        constraints.betas_GRB_tot.append(beta)
        constraints.Omega_GRB_tot.append(y)
    
    betas_grb1 = np.array(betas_grb1)
    M_grb1 = np.array(M_grb1)
    betas_grb2 = np.array(betas_grb2)
    M_grb2 = np.array(M_grb2)
    M_grb1_bbn = np.array(M_grb1_bbn)
    Omegas_grb1 = np.array(Omegas_grb1)
    M_grb2_bbn = np.array(M_grb2_bbn)
    Omegas_grb2 = np.array(Omegas_grb2)

    
    return M_grb1, M_grb2, betas_grb1, betas_grb2, Omegas_grb1, Omegas_grb2

def Betas_Reio(M_tot, omega):
    """
    Calculates the abundance of PBHs with a lifetime greater than the age of the universe leaves an imprint on the CMB through modifications of the ionization history and the damping of CMB anisotropies.

    Parameters:
            - M_tot (array-like): Array of masses in grams.
            - omega (float): This value is to assign the equation of state.

    Returns:
        A tuple containing three numpy arrays:
            -  M_reio (numpy.ndarray): Represents the masses of PBHs in the interval ( 1e15 , 1e17) grams.
            -  betas_reio (numpy.ndarray): Corresponds to the abundance obtained from M_reio.
            -  Omegas_reio (numpy.ndarray): Corresponds to evolution of abundance of betas_reio.
    """
    betas_reio = []
    M_reio = []

    M_reio_bbn = []
    Omegas_reio = []
    
    rho_form_rad = rho_f(M_tot, omega)

    for i in range(len(M_tot)):
        if M_tot[i]> 10**15 and M_tot[i]<10**17:
            M_reio.append(M_tot[i])
            beta = 2.4*10**(-26)*(M_tot[i]/(4.1*10**14))**4.3/constants.gam_rad**(1/2)
            betas_reio.append(beta)
        
            ln_den_f = np.log(rho_form_rad[i])
            if ln_den_f <= ln_den_end:
                continue
            ln_den = np.linspace(ln_den_f,ln_den_end,10000)
            sol_try = solve_ivp(diff_rad,(ln_den_f,ln_den_end),np.array([1.,0.]),events=end_evol,t_eval=ln_den,args=(M_tot[i],beta),method = "DOP853")
            Delta_t = constants.t_pl*(M_tot[i]/constants.M_pl_g)**3
            y = beta*sol_try.y[0][-1]*(1.-sol_try.y[1][-1]/Delta_t)**(1./3)
            Omegas_reio.append(y)
            M_reio_bbn.append(M_tot[i])
        else:
            beta = constants.ev1
            y = constants.ev2
        constraints.betas_Reio_tot.append(beta)
        constraints.Omega_Reio_tot.append(y)
    
    betas_reio = np.array(betas_reio)
    M_reio = np.array(M_reio)
    M_reio_bbn = np.array(M_reio_bbn)
    Omegas_reio = np.array(Omegas_reio)
    
    return M_reio, betas_reio, Omegas_reio

#def Betas_LSP(M_tot):

#    betas_lsp = []
#    M_lsp = []
    
#    for i in range(len(M_tot)):
#        if M_tot[i]< 10**11:
#            M_lsp.append(M_tot[i])
#            beta = 10**(-18)*(M_tot[i]/(10**11))**(-1/2)/constants.gam_rad**(1/2)
#            betas_lsp.append(beta)
#        else:
#            beta = constants.ev1
#        constraints.betas_LSP_tot.append(beta)
    
#    betas_lsp = np.array(betas_lsp)
#    M_lsp = np.array(M_lsp)
    
#   return M_lsp, betas_lsp


def Betas_LSP(M_tot, w):
    """
    Calculates the abundance of PBHs that may produce the lightest supersymmetric particles (LSP), predicted in supersymmetry and supergravity models, which are stable and may contribute to the totality of the DM in the universe.

        Parameters:
            - M_tot (array-like): Array of masses in grams.
            - w (omega) (float): This value is to assign the equation of state.

    Returns:
        A tuple containing three numpy arrays:
            -  M_lsp (numpy.ndarray): Represents the masses of PBHs < 1e11( m_lsp / 100 GeV)^(-1) grams.
            -  betas_lsp (numpy.ndarray): Corresponds to the abundance obtained from M_lsp.
            -  Omegas_lsp_tot (numpy.ndarray): Corresponds to evolution of abundance of betas_lsp.
    """

    ev1 = 1e-5
    ev2 = 1e-2
    M_pl = 1.22089*10**19

    t_pl_s = 5.39 * 10**-44
    s_to_evm1 = (1./6.5823) * 10**25
    t_pl = t_pl_s * s_to_evm1

    gam_rad = (1./3)**(3./2)
    H_end = 4.44 * 10**13.
    rho_end_inf = 3. * M_pl**2. * H_end**2.

    betas_lsp = []
    betas_lsp_tot = []
    M_lsp = []
    M_lsp_bbn = []
    Omegas_lsp = []
    M_lsp_pbbn = []
    Omegas_lsp_pbbn = []
    Omegas_lsp_tot = []

    #k_end_over_k_rad = (M_tot * H_end / (gam_rad * 3 * M_pl))**(1 + 3 * w) / (3 * (1 + w))
    #rho_form_rad = rho_end_inf / (k_end_over_k_rad)**((6 * (1 + w)) / (1 + (3 * w)))
    rho_form_rad = rho_f(M_tot, w)
    rho_end = (1e-2)**4
    ln_den_end = np.log(rho_end)

    for i in range(len(M_tot)):
        if M_tot[i] < 10**11:
            M_lsp.append(M_tot[i])
            beta = 10**(-18) * (M_tot[i] / (10**11))**(-1/2) / gam_rad**(1/2)
            betas_lsp.append(beta)

            ln_den_f = np.log(rho_form_rad[i])
            if ln_den_f <= ln_den_end:
                continue
            ln_den = np.linspace(ln_den_f, ln_den_end, 10000)
            sol_try = solve_ivp(diff_rad, (ln_den_f, ln_den_end), np.array([1., 0.]), events=end_evol, t_eval=ln_den, args=(M_tot[i], beta), method="DOP853")
            if sol_try.t[-1] > ln_den_end:
                sol_try = solve_ivp(diff_rad_rel, (ln_den_f, ln_den_end), np.array([1.]), t_eval=ln_den, args=(M_tot[i], beta), method="DOP853")
                y = beta * sol_try.y[0][-1] * (constants.M_pl_g / M_tot[i])
                Omegas_lsp_pbbn.append(y)
                M_lsp_pbbn.append(M_tot[i])
            else:
                Delta_t = t_pl * (M_tot[i] / constants.M_pl_g)**3
                y = beta * sol_try.y[0][-1] * (1. - sol_try.y[1][-1] / Delta_t)**(1./3)
                Omegas_lsp.append(y)
                M_lsp_bbn.append(M_tot[i])
        else:
            beta = ev1
            y = ev2
        constraints.betas_LSP_tot.append(beta)
        constraints.Omega_LSP_tot.append(y)
        betas_lsp_tot.append(beta)
        Omegas_lsp_tot.append(y)

    betas_lsp = np.array(betas_lsp)
    betas_lsp_tot = np.array(betas_lsp_tot)
    M_lsp = np.array(M_lsp)
    M_lsp_bbn = np.array(M_lsp_bbn)
    Omegas_lsp = np.array(Omegas_lsp)
    M_lsp_pbbn = np.array(M_lsp_pbbn)
    Omegas_lsp_pbbn = np.array(Omegas_lsp_pbbn)

    return M_lsp, betas_lsp, Omegas_lsp_tot


def get_Betas_full(M_tot):
    """
    This function calculates composite constraint values derived from various PBH constraints.
        Parameters:
                - M_tot (array-like): Array of masses in grams.

        Returns:
                - betas_full (numpy.ndarray): Represent the most robust constraints across diverse scenarios for each specific mass value. This output is saved in the module called constraints into variable named betas_full.
    """

    DM_tot = np.array(constraints.betas_DM_tot)
    BBN_tot = np.array(constraints.betas_BBN_tot)
    SD_tot = np.array(constraints.betas_SD_tot)
    CMB_tot = np.array(constraints.betas_CMB_AN_tot)
    GRB_tot = np.array(constraints.betas_GRB_tot)
    Reio_tot = np.array(constraints.betas_Reio_tot)
    LSP_tot = np.array(constraints.betas_LSP_tot)

    constraints.betas_full = M_tot * 0

    for i in range(len(M_tot)):
        # Collect only the non-empty arrays
        values = []
        if DM_tot.size: values.append(DM_tot[i])
        if BBN_tot.size: values.append(BBN_tot[i])
        if SD_tot.size: values.append(SD_tot[i])
        if CMB_tot.size: values.append(CMB_tot[i])
        if GRB_tot.size: values.append(GRB_tot[i])
        if Reio_tot.size: values.append(Reio_tot[i])
        if LSP_tot.size: values.append(LSP_tot[i])

        if values:  # Ensure there are values to calculate the minimum
            constraints.betas_full[i] = min(values)

    return constraints.betas_full


def get_Omegas_full(M_tot):
    """
    This function calculates composite constraint values derived from various PBH constraints.
        Parameters:
                - M_tot (array-like): Array of masses in grams.

        Returns:
                - Omegas_full (numpy.ndarray): Represent the most robust constraints across diverse scenarios for each specific mass value. This output is saved in the module called constraints into variable named Omegas_full.
    """

    DM_tot = np.array(constraints.Omega_DM_tot)
    BBN_tot = np.array(constraints.Omega_BBN_tot)
    SD_tot = np.array(constraints.Omega_SD_tot)
    CMB_tot = np.array(constraints.Omega_CMB_AN_tot)
    GRB_tot = np.array(constraints.Omega_GRB_tot)
    Reio_tot = np.array(constraints.Omega_Reio_tot)
    LSP_tot = np.array(constraints.Omega_LSP_tot)

    constraints.Omegas_full = M_tot * 0

    for i in range(len(M_tot)):
        # Collect only the non-empty arrays
        values = []
        if DM_tot.size: values.append(DM_tot[i])
        if BBN_tot.size: values.append(BBN_tot[i])
        if SD_tot.size: values.append(SD_tot[i])
        if CMB_tot.size: values.append(CMB_tot[i])
        if GRB_tot.size: values.append(GRB_tot[i])
        if Reio_tot.size: values.append(Reio_tot[i])
        if LSP_tot.size: values.append(LSP_tot[i])

        if values:  # Ensure there are values to calculate the minimum
            constraints.Omegas_full[i] = min(values)

    return constraints.Omegas_full



def inverse_error(betas, delta_c):
    aux = []
    for i in range(len(betas)):
        aux.append(delta_c/(np.sqrt(2)*special.erfcinv(betas[i])))
    return aux


def a_endre(rho_r0, rho_end_re):
    return (rho_r0 / rho_end_re) ** (1. / 4)


def k_rad(M):
    a_end_inf_rad = (constants.rho_r0 / constants.rho_end_inf) ** (1. / 4)
    k_end = a_end_inf_rad * constants.H_end
    k_end_over_k_rad = (M/(7.1*10**-2*constants.gam_rad*(1.8*10**15/constants.H_end)))**(1/2)
    k = (k_end/k_end_over_k_rad)*constants.GeV*constants.metter_m1
    k = np.array(k)
    return k