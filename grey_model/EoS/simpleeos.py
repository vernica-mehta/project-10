''''
This is me trying to walk through how eos.py works for a simple case of Hydrogen 
To be used for presentation since I have had trouble if the Regular EOS to get it to work with molecules. 

'''

import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u
import scipy.optimize as op
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import astropy.io.fits as pyfits
from astropy.table import Table


# Reading in Equilibrium coeffiencents from Tsuji_K
tsuji_K = Table.read('/Users/griffinkatrivesisbrown/Library/Mobile Documents/iCloud~md~obsidian/Documents/project-10/grey_model/data/tsuji_K.dat', format='ascii.fixed_width_no_header', \
 names=('mol', 'c0', 'c1', 'c2','c3','c4','molcode','ediss','comment'), col_starts=(0,7,19,31,43,55,67,80,87))

# Constants from eos.py 
debroglie_const = (c.h**2/2/np.pi/c.m_e/c.k_B).cgs.value
eV_Boltzmann_const = (u.eV/c.k_B).cgs.value
deybe_const = (c.k_B/8/np.pi/c.e.esu**2).cgs.value
delchi_const = (c.e.esu**2/(1*u.eV)).cgs.value

# I simplified the compostion function to work with just H 

def composition():
    """Return hydrogen properties only"""
    elt_names = np.array(['H'])
    n_p = np.array([1])
    masses = np.array([1.0])
    abund = np.array([1.0])  # 100% H
    ionI = np.array([13.595])
    ionII = np.array([-0.754])  # H- (negative ion)
    gI = np.array([2])
    gII = np.array([1])
    gIII = np.array([1])
    return abund, masses, n_p, ionI, ionII, gI, gII, gIII, elt_names

# Now we need to create a saha equation first as this will give us an idea of how many elements are ionised 

"""

Currently Ignoring all Molecules expect H2 

"""


def saha_h2_only(n_e, T):
    """
    Full Saha equation solver for H only taken from eos.saha,
    
    Solves matrix system for just H, H+, H- number densities
    Includes Debye screening and proper matrix formulation
    
    Debye screening -> Tells us how electric fields weakened in a plasma or ionized gas
    due to the collective response of free charged particles.

        This enables the problem to be a simple Ax=b linear problem.
    Results from this function can be used to solve the Saha equation as e.g. a function 
    of rho and T via e.g. tabulating or solving.

    Parameters
    ----------
    n_e: double
        the electron number density in cm^{-3}
    T: double
        Temperature in K.
    
    Returns
    -------
    rho: astropy.units quantity compatible with density
    mu: Mean molecular weight (dimensionless, i.e. to be multiplied by the AMU)
    ns: A vector of number densities of H, H+, H-, in cm^{-3}
    
    """
    
    # loads in our composition 
    abund, masses, n_p, ionI, ionII, gI, gII, gIII, elt_names = composition()
    n_elt = len(n_p)  # Just 1 for H
   
    
    # Debye screening 
    #according to Mihalas 9-178
    
    deybe_length = np.sqrt(deybe_const*T/n_e)
    z1_delchi = delchi_const/deybe_length
    
    # Low temperature case 
    # Ionisation fraction to zero 
    if T < 2000:
        ns = np.zeros(n_elt * 3)
        ns[3*np.arange(n_elt)] = abund  # Neutral atoms
        ns = np.maximum(n_e * 1e15 * ns, 1e-300)
    else:
        # Thermal de Broglie wavelength 
        debroglie = np.sqrt(debroglie_const / T)
        
        # Construct matrix 
        A = np.zeros((3*n_elt, 3*n_elt));
        for i in range(n_elt):
            # First equation: abundance ratio 
            # Sum of all H states / Sum of all H states = abund[i]/abund[0]
            #densities to H number densities is the ratio of abundances.
            A[3*i, 3*i:3*(i+1)] = [-abund[0], -abund[0], -abund[0]]
            A[3*i, :3] = [abund[i], abund[i], abund[i]]
           

            # Single ionization equilibrium for element
            he1 = 2./debroglie**3 *gII[i]/gI[i]*np.exp(-(ionI[i] - n_p[i]*z1_delchi)*eV_Boltzmann_const/T)
            A[3*i+1, 3*i:3*i+2] = [-he1/n_e, 1]  # -K/n_e * n_neutral + n_ion = 0
            # H- formation (special case for i=0)
            if i == 0:
                he2 = 2. / debroglie**3 * gI[i] / gIII[i] * np.exp(-(np.abs(ionII[i]) - n_p[i]*z1_delchi) * eV_Boltzmann_const / T)
                A[3*i+2, 3*i:3*i+3] = [1, 0, -he2/n_e]  # n_neutral - K/n_e * n_Hminus = 0
            else:
                # Double ionization (not relevant for H only, but keeping structure)
                he2 = 2. / debroglie**3 * gIII[i] / gII[i] * np.exp(-(ionII[i] - n_p[i]*z1_delchi) * eV_Boltzmann_const / T)
                A[3*i+2, 3*i+1:3*i+3] = [-he2/n_e, 1]
            

        #Add in the equation for electron number density. This over-writes the equation
        #n_H + n_H  
        A[0, :] = np.concatenate(([0, 1, -1], np.tile([0, 1, 2], n_elt-1)))
        
        #where b only has one non-zero number (the electron 
        #number density)
        b = np.zeros(3*n_elt)
        b[0] = n_e
        ns = np.linalg.solve(A, b)
        
        # Make sure all densities are postive  
        ns = np.maximum(ns, 1e-300)
    
    # High temperature/pressure ionization correction 
    ns_highT = np.zeros(n_elt * 3)
    ns_highT[1] = abund[0]  # Fully ionized H
    ns_highT[2 + np.arange(n_elt-1)*3] = abund[1:]  # Other elements ionized
    ns_highT = ns_highT / (abund[0] + 2*np.sum(abund[1:])) * n_e
    atom_size = 1e-8  # cm
    if n_e * atom_size**3 > 2:
        ns = ns_highT
    elif n_e * atom_size**3 > 1:
        frac = (n_e * atom_size**3 - 1) / 1.0
        ns = frac * ns_highT + (1 - frac) * ns
    
    # Calculate total H nuclei 
    n_h = np.sum(ns[:3])
    
    # Density 
    rho_in_g_cm3 = n_h * np.sum(abund * masses) * c.u.to(u.g).value
    
    #Fractional "abundance" of electrons.
    f_e = n_e / n_h
    # Mean molecular weight 
    mu = np.sum(abund * masses) / (np.sum(abund) + f_e)
    
    # Specific internal energy: ionization energy per unit mass
    # Normalized per H nucleus (reference unit) but accounting for 
    #  the mean molecular weight of the full gas mixture 
    Ui = (ns[1] * 13.6) * u.eV / n_h / np.sum(abund * masses * u.u)
    
    return rho_in_g_cm3, mu, Ui, ns



def saha_solve_h2(log_n_e_mol_cm3, T_K, rho_0_in_g_cm3):
    """
    Find the value of n_e (electron density) that produces the correct total mass density
    (taken from eos.saha_solve)
    """
    n_e = np.exp(log_n_e_mol_cm3[0]) * c.N_A.value
    rho, mu, Ui, ns = saha_h2_only(n_e, T_K)
    return np.log(rho_0_in_g_cm3 / rho)

def ns_from_rho_T_h2(rho, T):
    """
    Compute ionization state from density and temperature (like eos.ns_from_rho_T)
    
    Uses full Saha equation to get H, H+, H- abundances ignoring molecules.
    """
    rho_in_g_cm3 = rho.to(u.g/u.cm**3).value

    
    # Initial guess for electron density
    x0 = np.log(rho_in_g_cm3)
    T_K = np.maximum(T.to(u.K).value,1000)
    x0 += np.log(2/(10*np.exp(40e3/T_K) + 1))
    
    # Solve for n_e that gives correct density 
    print(T.cgs.value)
    res = op.fsolve(saha_solve_h2, x0, args=(T.cgs.value, rho_in_g_cm3), xtol=1e-6)
    n_e = np.exp(res[0]) * c.N_A.value
    
    
    # Get final results from full Saha equation
   
    rho_check, mu, Ui, ns = saha_h2_only(n_e, T.cgs.value)
    
    # Verify convergence
    print(T_K)
    print(T.cgs.value) 
    if np.abs(rho_check / rho_in_g_cm3 - 1) > 0.01:
        print(f"Warning: Density check failed ({rho_check:.3e} vs {rho_in_g_cm3:.3e})")
    
  
    
    return n_e*(u.cm**(-3)), ns*(u.cm**(-3)), mu, Ui


"""
Dealing with molecules 

"""

def equilibrium_equation_h2_only(rho, T):
    """
    Simplified version of eos.equilibrium_equation() for H2 only.
    """
    abund, masses, n_p, ionI, ionII, gI, gII, gIII, elt_names = composition()
    
    # Reference pressure,which is the ideal gas pressure corresponding to
    #atomic Hydrogen only. 
    log_P_ref = np.log10((rho * c.k_B * T / u.u).to(u.dyn/u.cm**2).value)
    
    # Get H2 data
    nmol = 1 # we only have one molecule 
    natom = len(abund)

   
    '''
    In original it was: 
    The linear matrix, with one equation (row) per atom,  one
    for the electron.
    Species ordering: [e-, neutrals (natom), ions (natom), molecules (nmol)]
     
    For us we now only need to worry about H so: 
    '''
    # For H only: [e-, H, H+, H2]
    # Total columns = 1 + 2*natom + nmol = 1 + 2*1 + 1 = 4
    # Linear matrix (conservation equations)
    # Rows: 1 (electron) + natom (atom conservation) = 1 + 1 = 2 
    linear_matrix = np.zeros((natom+1, 1 + 2*natom + nmol))
    linear_b = np.zeros((natom + 1)) 


    # Logarithmic matrix (equilibrium equations) 
    # Rows: natom (ionization) + nmol (molecules) = 1 + 1 = 2
    log_matrix = np.zeros((natom + nmol, 1 + 2*natom + nmol))
    log_b = np.zeros(natom + nmol)

    #The theta value for computing the molecular equilibrium constants
    theta = 5040/T.to(u.K).value
    #The equivalent for the  atoms
    eV_kTln10 = float(1*u.eV/c.k_B/T/np.log(10))

    #First, the electron equation
    # All ions contribute: -P_e + sum(P_ions) = 0
    linear_matrix[0,1+natom:1+2*natom] = 1 # Original would be all ions but for single H is just H+ 
    linear_matrix[0,0] = -1 # Negative electron pressure 

    #Next, the RHS of the equation
    linear_b[1:] = abund/np.sum(masses*abund)

    #Prod in = Prod loss 
    # Total H atoms (in) = Total H atoms (out)
    # Now, the atomic part of the equation 
    for i in range(natom):
        linear_matrix[i+1, i+1] = 1           # Neutral atom
        linear_matrix[i+1, natom+i+1] = 1     # Ion

    

    debroglie = np.sqrt(debroglie_const / T.to(u.K).value)
    kBT = (c.k_B * T).cgs.value
    
    #The logarithmic part of the matrix for atoms.
    for i in range(natom):
        log_b[i] = np.log10(2*kBT/debroglie**3 *gII[i]/gI[i]) - eV_kTln10*ionI[i]
        log_matrix[i,0] = 1
        log_matrix[i,1+i] = -1
        log_matrix[i,1+i+natom] = 1

   
    #Next, the molecular part of the equation
    for i in range(nmol):
        log_matrix[natom + i, 2*natom+1+i] = -1 # -log(P_molecule)

        # This is not very useful in this case since we are only dealing with one molecule, mol = HH
        # If there was more molecules it would go through each one
        mol = tsuji_K[i]
        
        # Applying chemical equilibrium coeffiencents from Tsuji_K
        log_b[natom + i] = mol['c0'] + mol['c1']*theta + mol['c2']*theta**2 + mol['c3']*theta**3 + mol['c4']*theta**4

        molcode = mol['molcode']
        n_element_types = int(molcode[0]) # Number of different elements in molecule

        for j in range(n_element_types):
            atom = int(molcode[1+j*3:3+j*3])
            natom_in_mol = int(molcode[3+j*3:4+j*3])
            k = np.argwhere(n_p==atom)[0,0]
            linear_matrix[k+1,2*natom+1+i] = natom_in_mol
            log_matrix[natom + i, 1+k] = natom_in_mol
    
    return linear_matrix, linear_b, log_matrix, log_b, log_P_ref



def eq_solve_func_h2(logps, linear_matrix, linear_b, log_matrix, log_b, log_P_ref, abund):
    """
    Simplified version of eos.eq_solve_func() 
    
    Combines logarithmic equilibrium equations with linear conservation equations.

    Returns:
        --------
        resid: np.array
            concatenated log resids (n_atoms + n_mol) then concatenated linear resids 
            (1 + n_atoms)
    """
    ps = 10**(logps - log_P_ref)
    
    # Logarithmic part (equilibrium)
    log_part = np.dot(log_matrix, logps) - log_b
    
    # Linear part (conservation)
    linear_part = np.dot(linear_matrix, ps) - linear_b
    #Dividing the linear part by the abundance just makes all numbers near 0 and the 
    #same order of magnitude.
    linear_part[1:] /= abund  # Normalize
    
    return np.concatenate((log_part, linear_part))


def equilibrium_solve_h2(rho, T):
    """
    Simplified version of eos.equilibrium_solve() 
    
    Solves for equilibrium of H, H+, e-, and H2 at given density and temperature.
    
    Returns:
    --------
    log(pressure) in dyn/cm^2 for: [e-, H, H+, H2]
    """
    linear_matrix, linear_b, log_matrix, log_b, log_P_ref = equilibrium_equation_h2_only(rho, T)
    abund, masses, n_p, ionI, ionII, gI, gII, gIII, elt_names = composition()
    
    
    #Starting point -> using saha equation to get our initial guess 
    
    n_e, ns, mu, Ui = ns_from_rho_T_h2(rho,T)
    
   
     #Convert number density to log partial pressure.
    if min(ns)<0:
        import pdb; pdb.set_trace()
    logp_atom = np.log10((ns*c.k_B*T).to(u.dyn/u.cm**2).value)
    logpe = np.log10((n_e*c.k_B*T).to(u.dyn/u.cm**2).value)    
    x0 = np.ones(linear_matrix.shape[1])*log_P_ref - 12
   
   
    x0[0] = logpe #n_e
    x0[1+np.arange(len(abund))] = logp_atom[3*np.arange(len(abund))] #neutrals
    x0[len(abund)+1+np.arange(len(abund))] = logp_atom[3*np.arange(len(abund))+1] #ions

    # Initial guess 
    #x0 = np.ones(linear_matrix.shape[1])*log_P_ref - 12
    #x0[1:len(abund)+1] = log_P_ref + np.log10(abund)

    
    # Start with mostly neutral H
    x0[1] = log_P_ref  # H (neutral)
    x0[0] = log_P_ref - 6  # e- (low)
    x0[2] = log_P_ref - 6  # H+ (low)
    
    # H2 depends on temperature
    if T.to(u.K).value < 2000:
        x0[3] = log_P_ref - 0.5  # Mostly H2 at low T
        x0[1] -= 2  # Less atomic H
    else:
        x0[3] = log_P_ref - 5  # Little H2 at high T
    
    
    res = op.root(eq_solve_func_h2, x0, 
                  args=(linear_matrix, linear_b, log_matrix, log_b, log_P_ref, abund),
                  method='lm')
    
    if res.success:
        return res.x
    else:
        print("WARNING: Solution did not converge")
        return res.x


def eos_rho_T_h2_only(rho, T):
    """
    Big approximation because dealing the saha equation is a pain:
    doesn't include gamma, mu , Ui or ne 
    
    Returns:
    --------
    P: Total pressure
    species_pressures: Dictionary with partial pressures
    """
    # Solve equilibrium
    logps = equilibrium_solve_h2(rho, T)
    
    # Convert to pressures
    ps = 10**logps * u.dyn / u.cm**2
    
    species_pressures = {
        'e-': ps[0],
        'H': ps[1],
        'H+': ps[2],
        'H2': ps[3]
    }
    
    # Total pressure
    P_total = np.sum(ps)
    
    return P_total, species_pressures




if __name__ == '__main__':

    print("SIMPLIFIED EOS.PY - H2 ONLY VERSION")
 
    
    # Test at different temperatures
    rho = 1e-7 * u.g / u.cm**3
    temperatures =  [100,500,1000,2000,3000,5000, 8000] * u.K
    
    print(f"\nAt density rho = {rho}:")
    print(f"{'T (K)':<10} {'P_total':<15} {'P(H)':<15} {'P(H2)':<15} {'Dominant'}")
    print("-"*70)
    
    for T in temperatures:
        P_total, species = eos_rho_T_h2_only(rho, T)
        print(P_total)
            
        dominant = 'H2' if species['H2'] > species['H'] else 'H'
            
        print(f"{T.value:<10.0f} {P_total.value:<15.2e} "
                f"{species['H'].value:<15.2e} {species['H2'].value:<15.2e} {dominant}")
    
    # Make a comparison plot
    print("\nGenerating comparison plot...")
    T_range = np.logspace(2.8, 4.0, 30) * u.K
    
    P_H_list = []
    P_H2_list = []
    P_total_list = []
    
    for T in T_range:
        P_total, species = eos_rho_T_h2_only(rho, T)
        P_H_list.append(species['H'].value)
        P_H2_list.append(species['H2'].value)
        P_total_list.append(P_total.value)

    plt.figure(figsize=(12, 5))
    
    # Plotting Pressure vs Temp
    plt.subplot(1, 2, 1)
    plt.loglog(T_range.value, P_H_list, 'r-', linewidth=2, label='P(H)')
    plt.loglog(T_range.value, P_H2_list, 'b-', linewidth=2, label='P(H2)')
    plt.loglog(T_range.value, P_total_list, '--', linewidth=1.5, label='P(total)')
    plt.xlabel('Temperature (K)')
    plt.gca().set_xscale('linear')
    plt.ylabel('Pressure (dyne/cm2)')
    plt.title('Partial Pressures, no gamma, mu , Ui or ne calculated ')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plotting how many H we get:
    plt.subplot(1, 2, 2)
    fraction_H2 = np.array(P_H_list) / (np.array(P_H_list) + np.array(P_H2_list))
    plt.semilogx(T_range.value, fraction_H2, 'g-', linewidth=2)
    plt.gca().set_xscale('linear')
    plt.xlabel('Temperature (K)')
    plt.ylabel('H/H_total')
    plt.title('Dissociation Fraction')
    plt.grid(True, alpha=0.3)
    plt.ylim([0, 1])
    
    plt.tight_layout()
    plt.savefig('h2_eos_simple.png', dpi=150, bbox_inches='tight')
    print("Plot saved as 'h2_eos_simple.png'")
    plt.show()