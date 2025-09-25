"""
Non-Grey Atmospheric Model
Based on the work of Parmentier & Guillot (2014)

Model will produce a temperature-optical depth profile for a given set of arguments,
arguments can be changed to see how the temperature profile changes.

"""

# Import visible opacities from grey_model.py by Vernica to use in non-grey model
# -- Import Packages -- #
import numpy as np
import astropy.constants as c
import astropy.units as u
import matplotlib.pyplot as plt
import sys
import argparse

from grey_model import GreyModel
model = GreyModel()

model.apply_opacs()

kappa_ross = model.kappa_ross
kappa_v = model.kappa_visual

gamma_v = kappa_v / kappa_ross





"""
Non-Grey Atmospheric Model
Based on the work of Parmentier & Guillot (2014)

Model will produce a temperature-optical depth profile for a given set of arguments,
arguments can be changed to see how the temperature profile changes.

"""

# Import visible opacities from grey_model.py by Vernica to use in non-grey model
from grey_model import GreyModel
model = GreyModel()
ross_opac = model.kappa_ross
visual_opac = model.kappa_visual

print(ross_opac)





# -- File Arguments -- #

parser = argparse.ArgumentParser()  # Sets up argparse package to define and check for arguments

parser.add_argument("--T_eff", required=False, type=int, default=143) # Effective surface temperature in Kelvin  
parser.add_argument("--T_eff_sun", required=False, type=int, default=6500) # Irradiation temperature in Kelvin
parser.add_argument("--mu", required=False, type=float, default=1) # Cosine of angle of incident radiation, default is 1 (normal incidence)

parser.add_argument("--modelParamChange", required=False, type=bool, default=True) # If True: Changing thermal opacities, variables become R and beta.
                                                                                   # If False: Changing overall profile, variables become kappa_P and tau_lim
                                                                                   
parser.add_argument("--R", required=False, type=float, default=0.1) # Ratio of Thermal line and continuum opacities kappa_1/kappa_2
parser.add_argument("--beta", required=False, type=float, default=0.5) # Partitioning of thermal opacity between two channels 

parser.add_argument("--tau_lim", required=False, type=float, default=100) # Optical depth limit for picket fence model
parser.add_argument("--gamma_P", required=False, type=float, default=0.4) # Ratio of Planck mean thermal opacity to Roeseland mean visible opacity kappa_P/kappa_v

args = parser.parse_args() # Checks if set up args are present in run line, if not sets default values




# -- Call Arguments -- #
kwargs = {
    "T_eff": args.T_eff*u.K,
    "T_eff_sun": args.T_eff_sun*u.K,
    "mu": args.mu,
    "gamma_v":gamma_v,
    "modelParamChange": args.modelParamChange,
    "R": args.R,
    "beta": args.beta,
    "tau_lim": args.tau_lim,
    "gamma_P": args.gamma_P,
    "d": args.d,
    "R_sun": args.R_sun
}





# -- Constants -- #
h = c.h.cgs.value # Planck's constant in cgs
k_B = c.k_B.cgs.value # Boltzmann's constant in cgs
c = c.c.cgs.value # Speed of light in cgs
sigma_SB = c.sigma_sb.cgs.value # Stefan-Boltzmann constant in cgs


def nonGreyModel(**kwargs):
    """
    kwargs
    -----------
    T_eff: float
        Effective temperature of planet in Kelvin.
        default is 144K (Jupiter-like)
    T_eff_sun : float
        Effective temperature of sun in Kelvin.
        default is 6500K (solar irradiation)
    kappa_v : array
        Visible opacity in cm^2/g

  
    modelParamChange : bool
        if True: Changing thermal opacities, variables become R and beta
        if False: Changing overall profile, variables become kappa_P and tau_lim
    
    R : array
        Ratio of Thermal line and continuum opacities kappa_1/kappa_2
        Only used if modelParamChange is True
    beta : array
        Partitioning of thermal opacity between two channels
        Only used if modelParamChange is True

        
    tau_lim : float
        Optical depth limit for picket fence model
        Only used if modelParamChange is False
    gamma_P : float
        Ratio of Planck mean thermal opacity to Roeseland mean visible opacity kappa_P/kappa_v
        Only used if modelParamChange is False
    """
    # ----------------------------------------------------------------

    # -- Parameters -- #
    # Call constants from argumets
    gamma_v = kwargs.get('gamma_v') # Visible opacity
    

    # Based on which paramters are being changed for the model, define the other parameters accordingly
    if kwargs.get('modelParamChange', True) == True:
    
        # Get parameters from kwargs or set defaults
        
        beta = kwargs.get('beta', 0.5)
        R = kwargs.get('R', 0.1)

        # Use defined params to calculate other params for model
        planck_function = (2*h*nu**3/c.c**2)/(np.exp(h*nu/(k_B*T))-1)

        gamma_P = beta + R - beta*R + (beta + R - beta*R)/R - (beta + R - beta*R)**2/R
        tau_lim = (np.sqrt(R)*np.sqrt(beta*(R-1)**2 - beta**2 *(R - 1)**2 + R))/(np.sqrt(3)*(beta + R - beta*R)**2)
    
    else:
        # Get parameters from kwargs or set defaults
        gamma_P = kwargs.get('gamma_P', 0.4)
        tau_lim = kwargs.get('tau_lim', 100)

        # Calculate other params for model
        delta = 3*gamma_P + 3*np.sqrt(gamma_P) * tau_lim*(2*np.sqrt(3)*gamma_P + 3*gamma_P**(3/2)*tau_lim - 4*np.sqrt(3))

        R = (np.sqrt(3*gamma_P) + 3*gamma_P*tau_lim + np.sqrt(delta)) / (np.sqrt(3*gamma_P) +3*gamma_P*tau_lim - np.sqrt(delta))
        beta = (np.sqrt(delta) - np.sqrt(3*gamma_P) + 3*gamma_P*tau_lim)/(2*np.sqrt(delta))
    
    # Calculate gamma_1 and gamma_2
    gamma_1 = beta + R - beta*R
    gamma_2 = (beta + R - beta*R)/R

    # -- Coeffients -- #
    # Calculate temperature profile coefficients
    T_eff_sun = kwargs.get('T_eff_sun') # Irradiation temperature in Kelvin
    W = (R/D)**2 * mu**2
    T_irr = W * T_eff_sun

    # Internal temperature in Kelvin


    # Calculate sub-coefficients for temperature profile coefficients (from page 8 and pretty ugly, sorry)
    A_t_1 = gamma_1**2 * np.log(1 + 1/(tau_lim*gamma_1))
    A_t_2 = gamma_2**2 * np.log(1 + 1/(tau_lim*gamma_2))
    
    a_0 = 1/(gamma_1) + 1/(gamma_2)

    a_1 = -1/(3*tau_lim**2) * ( (gamma_P)/(1-gamma_P) * (gamma_1 + gamma_2 - 2)/(gamma_1 + gamma_2) +
             (gamma_1 + gamma_2)*tau_lim - (A_t_1 + A_t_2)*tau_lim**2)
    
    b_0 = ((gamma_1*gamma_2)/(gamma_1 - gamma_2) * (A_t_1 - A_t_2)/3 - (gamma_1*gamma_2)**2/np.sqrt(3*gamma_P) 
            - (gamma_1*gamma_2)**3/((1-gamma_1) * (1-gamma_2) * (gamma_1 + gamma_2))   )**-1
    
    A = 1/3 * (a_0 + a_1*b_0)

    B = -1/3 * (gamma_1*gamma_2)**2*b_0/gamma_P


    for i in range(0,len(gamma_v[0,:])):
    
        A_v_1 = gamma_1**2 * np.log(1 + gamma_v/gamma_1)
        A_v_2 = gamma_2**2 * np.log(1 + gamma_v/gamma_2)

        a_2 = tau_lim**2/(gamma_P*gamma_v**2) * ( (3*gamma_1**2 - gamma_v**2) * (gamma_1 + gamma_2) -
            3*gamma_v*(6*gamma_1**2*gamma_2**2 - gamma_v**2*(gamma_1**2 + gamma_2**2)))/(1 - gamma_v**2*tau_lim)

        a_3 = -(tau_lim**2 * (3*gamma_1**2 - gamma_v**2) * (A_v_2 + A_v_1))/(gamma_P*gamma_v**3 * (1 - gamma_v**2*tau_lim**2))
        

        b_1 = gamma_1*gamma_2*(3*gamma_1**2 - gamma_v**2)*(3*gamma_2**2 
            - gamma_v**2)*tau_lim/(gamma_P*gamma_v**2*(gamma_v**2*tau_lim**2 - 1))
    
        b_2 = 3*(gamma_1 + gamma_2)*gamma_v**3/(gamma_P*gamma_v**2*(gamma_v**2*tau_lim**2 -1))

        b_3 = (A_v_2 - A_v_1)/(gamma_v*(gamma_1 - gamma_2))
    

        
    
        C = -1/3 * (b_0*b_1*(1 + b_2 + b_3)*a_1 + a_2 + a_3)

        D = 1/3 * (gamma_1*gamma_2)**2/gamma_P * b_0*b_1*(1 + b_2 + b_3)

        E = (3 - (gamma_v/gamma_1)**2)*(3 - (gamma_v/gamma_2)**2)/(9*gamma_v*((gamma_v*tau_lim)**2 - 1))

    

    # -- Temperature Profile -- #
    # Define optical depth array
    tau = np.logspace(1e-6, 1e6, 50)

    # Define temperature profile
    # NOTE: the irradiated part can be expanded for more optical opacities, sum that part for each gamma_v
    T_4 = 3*T_int**4/4 * (tau + A + B*np.exp(-tau/tau_lim)) + 3*T_irr**4* mu/4*(C + D*np.exp(-tau/tau_lim) + E*np.exp(-gamma_v*tau))
    T = (T_4)**(1/4)

    # -- Plotting -- #

    fig, ax = plt.subplots(figsize=(8,6) )