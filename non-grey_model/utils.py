# utils.py
# Utility functions for atmosphere modelling

# imports
import numpy as np
import astropy.constants as c
import astropy.units as u

def T_tau(tau_h, Teff, Tirr, kappa_j, kappa_b, f, r, D):
	"""Calculate temperature profile as a function of optical depth.
	
	PARAMETERS
	----------
	tau_h : `array`
		Array of optical depths.
	Teff : `float`
		Effective temperature of the planet in Kelvin.
	Tirr : `float`
		Irradiation temperature from the star in Kelvin.
	kappa_j : `float`
		J-band opacity in cm^2/g.
	kappa_b : `float`
		Visible-band opacity in cm^2/g.
	f : `float`
		Redistribution factor (1/2 for dayside average, 1/4 for full average).
	r : `float`
		Radius of the star in meters.
	D : `float`
		Distance from the star to the planet in meters.

	RETURNS
	-------
		Temperature profile in Kelvin.
	"""  

	W = f * (r/D)**2 # dilution factor
	brackets = tau_h + 1/np.sqrt(3) + 4*W/3 * (Tirr/Teff)**4

	if Teff == 0: # handle case of no internal heat
		T = (kappa_j/kappa_b) * W * Tirr**4 * np.ones_like(tau_h)
	else:
		T = (kappa_j/kappa_b) * 3/4 * Teff**4 * (brackets)

	return (T**0.25).cgs

def dP_dTau(tau, P, tau_grid, T_grid, g, f_kappa):
	"""Calculate dP/dTau for hydrostatic equilibrium.
	
	PARAMETERS
	----------
	tau : `float`
		Optical depth at the layer.
	P : `float`
		Pressure at the layer in dyn/cm^2.
	tau_grid : `array`
		Array of optical depths for interpolation.
	T_grid : `array`
		Array of temperatures for interpolation.
	g : `float`
		Surface gravity in cm/s^2.
	f_kappa : `callable`
		Function to compute opacity given log10(P) and T.
		
	RETURNS
	-------
		Pressure gradient with respect to optical depth in dyn/cm^2.
	"""

	# interpolate temperature and opacity
	T = np.interp(tau, tau_grid, T_grid)
	kappa_bar = f_kappa((np.log10(P), T))

	return g / kappa_bar
    
def planck(nu, T):
	"""Calculate the Planck function at frequency nu and temperature T.
	
	PARAMETERS
	----------
	nu : `Quantity`
		Frequency in Hz.
	T : `Quantity`
		Temperature in Kelvin.

	RETURNS
	-------
		Planck function value in erg/(cm^2 s Hz).
	"""

	nu = nu.to(u.Hz)
	p1 = (2*c.h*nu**3) / (c.c**2)
	p2 = np.expm1((c.h*nu)/(c.k_B*T))

	return (p1/p2)

def irradiation(r, D, nu, Tirr):
	"""Calculate the irradiation function at frequency nu and irradiation temperature Tirr.

	PARAMETERS
	----------
	r : `float`
		Radius of the star in meters.
	D : `float`
		Distance from the star to the planet in meters.
	nu : `Quantity`
		Frequency in Hz.
	Tirr : `Quantity`
		Irradiation temperature in Kelvin.

	RETURNS
	-------
		Irradiation function value in erg/cm^2.
	"""

	dist = (r/D)**2
	
	return dist * planck(nu, Tirr)