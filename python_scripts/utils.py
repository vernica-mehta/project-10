import numpy as np
import astropy.constants as c
import astropy.units as u

def T_tau(tau_h, Teff, Tirr, kappa_j, kappa_b, f, r, D):  

	W = f * (r/D)**2

	brackets = tau_h + 1/np.sqrt(3) + 4*W/3 * (Tirr/Teff)**4

	if Teff == 0:
		T = (kappa_j/kappa_b) * W * Tirr**4 * np.ones_like(tau_h)
	else:
		T = (kappa_j/kappa_b) * 3/4 * Teff**4 * (brackets)

	return (T**0.25).cgs

def dP_dTau(_, P, T, g, f_kappa):
	kappa_bar = f_kappa((np.log10(P), T))
	return g / kappa_bar
    
def planck(nu, T):

	p1 = (2*c.h*nu**3) / (c.c**2)
	p2 = np.exp((c.h*nu)/(c.k_B*T)) - 1

	return (p1/p2).cgs

def irradiation(r, D, nu, Tirr):
	dist = (r/D)**2
	return dist * planck(nu, Tirr).cgs