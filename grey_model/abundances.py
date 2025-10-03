import numpy as np
import astropy.units as u

def estimate_molecular_abundances(T, P):
    """
    Estimate molecular abundances for different atmospheric compositions
    
    Parameters:
    -----------
    T : float
        Temperature in K
    P : float
        Pressure in dyne/cm^2
    
    Returns:
    --------
    dict : Molecular number densities in cm^-3
    """
    
    # Total number density from ideal gas law
    n_total = P / (1.38e-16 * T)  # cm^-3
    
    abundances = {
        '1H2-16O': n_total * 1e-2,   # H2O mixing ratio ~10^-4
        '14N-1H3': n_total * 1e-2,   # NH3 mixing ratio ~10^-4  
        '12C-1H4': n_total * 1e-2,   # CH4 mixing ratio ~10^-6
    }

    print(f"[Abundances] T={T}K, P={P:.2e}, n_total={n_total:.2e} cm^-3")
    for mol, abund in abundances.items():
        mixing_ratio = abund / n_total
        print(f"[Abundances] {mol}: {abund:.2e} cm^-3 (mixing ratio: {mixing_ratio:.2e})")
    
    return abundances