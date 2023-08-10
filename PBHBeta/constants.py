# File: constants.py
# Module with physical constants

t_pl_s = 5.39*10**-44  # time in seconds (float)
s_to_evm1 = (1./6.5823)*10**25  # conversion factor from seconds to GeV^-1 (float)
t_pl = t_pl_s*s_to_evm1  # time in GeV^-1 (float)
M_pl = 1.22089*10**19  # Planck mass in GeV (float)
M_pl_g = 2.17645e-5  # Planck mass in grams (float)
grams_to_solar_mass = 5.02785e-34  # conversion factor from grams to solar masses (float)
gam_rad = (1./3)**(3./2)  # radiation constant gamma (float)
H_end = 4.44*10**13  # inflation energy in GeV (float)
rho_end = (1e-2)**4  # energy density at the end of inflation in GeV^4 (float)
rho_end_inf = 3.*M_pl**2.*H_end**2.  # energy density at the end of inflation in Planck units (float)
ev1 = 1e-5  # energy scale in GeV (float)
ev2 = 1e-2  # energy scale in GeV (float)
GeV = (1./1.97)*10**16 # m^-1
metter_m1 = (1./3.24078)*10**23 # Mpc^-1
rho_c = 8.07 * 10 ** -47  # GeV^4
Om_r0 = 8.4 * 10 ** -5.
rho_r0 = rho_c * Om_r0
A = (5+3/3)**2/(4*(1+1/3)**2)
B =  25/4
C = (8/2)**2/4
