# molecular_lines.py
import pandas as pd
import numpy as np
import astropy.constants as c
import astropy.units as u
from scipy.special import voigt_profile

class MolecularLines:
    def __init__(self, csv_files):
        """
        Initialize with list of CSV files containing molecular line data
        
        Parameters:
        -----------
        csv_files : list
            List of paths to CSV files with molecular line data
        """
        self.molecules = {}
        self.load_line_data(csv_files)
    
    def load_line_data(self, csv_files):
        """Load molecular line data from CSV files"""
        for csv_file in csv_files:
            # Extract molecule name from filename
            molecule_name = self.extract_molecule_name(csv_file)
            
            # Read the CSV file
            df = pd.read_csv(csv_file)
            
            # Store the data
            self.molecules[molecule_name] = {
                'frequencies': df['nu'].values * u.cm**-1,  # Convert to proper units
                'strengths': df['A'].values,  # Einstein A coefficients
                'line_strengths': df['S'].values,
                'lower_energy': df['E"'].values * u.cm**-1,
                'upper_g': df["g'"].values,
                'lower_g': df['g"'].values
            }
    
    def extract_molecule_name(self, filename):
        """Extract molecule name from filename"""
        # Based on your filenames like "20251002035654__14N-1H3__144.0K.csv"
        import os
        basename = os.path.basename(filename)
        parts = basename.split('__')
        if len(parts) >= 2:
            return parts[1]  # e.g., "14N-1H3" or "1H2-16O"
        return basename.split('.')[0]
    
    def compute_molecular_opacity(self, nu_grid, T, P, molecule_abundances, max_lines=100):
        """
        Compute molecular line opacity for given conditions
        
        Parameters:
        -----------
        nu_grid : array
            Frequency grid in Hz
        T : float
            Temperature in K
        P : float  
            Pressure in dyne/cm^2
        molecule_abundances : dict
            Dictionary of molecule abundances (number densities)
        """
        print(f"[MolecularLines] Computing opacity for T={T}K, P={P:.2e} dyne/cm^2")
        print(f"[MolecularLines] Frequency grid: {len(nu_grid)} points")

        total_opacity = np.zeros_like(nu_grid)
        
        for molecule, data in self.molecules.items():
            if molecule not in molecule_abundances:
                continue

            print(f"[MolecularLines] Processing {molecule}...")
                
            # Get line data
            line_freqs = (data['frequencies'] * c.c).to(u.Hz).value
            strengths = data['line_strengths']
            lower_energies = data['lower_energy']
            
            # Limit to strongest lines only
            n_total_lines = len(line_freqs)
            if n_total_lines > max_lines:
                sort_indices = np.argsort(strengths)[::-1][:max_lines]
                line_freqs = line_freqs[sort_indices]
                strengths = strengths[sort_indices]
                lower_energies = lower_energies[sort_indices]
            
            n_lines = len(line_freqs)
            print(f"[MolecularLines] {molecule}: Processing {n_lines} lines...")
            
            # Pre-calculate constants
            number_density_val = molecule_abundances[molecule].value if hasattr(molecule_abundances[molecule], 'value') else molecule_abundances[molecule]
            molecular_mass = self.get_molecular_mass(molecule)
            
            # Pre-calculate Boltzmann factors for all lines (vectorized)
            kT_cgs = c.k_B.cgs.value * T / u.K # erg
            E_lower_erg = lower_energies.to_value(u.cm**-1) * c.h.cgs.value * c.c.cgs.value  # erg
            boltzmann_factors = np.exp(-E_lower_erg / kT_cgs)
            
            for i in range(n_lines):
                nu0 = line_freqs[i]
                S = strengths[i]
                boltzmann = boltzmann_factors[i]
                
                # Line strength at temperature T
                line_strength_val = S * boltzmann
                if hasattr(line_strength_val, 'value'):
                    line_strength_val = line_strength_val.value
                
                # Doppler width
                molecular_mass = self.get_molecular_mass(molecule)
                doppler_width = (nu0 * u.Hz * np.sqrt(c.k_B * T / 
                                                    molecular_mass) / c.c).to(u.Hz).value
                
                # Broadening
                gamma_natural = 1e6
                gamma_pressure = P * 1e-3
                gamma_total = gamma_natural + gamma_pressure
                
                # Voigt profile
                x = (nu_grid - nu0) / doppler_width
                profile = voigt_profile(x, 1.0, gamma_total/doppler_width)
                profile = profile / (doppler_width * np.sqrt(2 * np.pi))  # Hz^-1
                
                # Calculate opacity contribution with proper physical units
                opacity_contribution = (number_density_val * line_strength_val * 
                                      profile * 3e10)
                
                total_opacity += opacity_contribution

            print(f"[MolecularLines] {molecule}: Done processing {n_lines} lines")

        return total_opacity
    
    def get_molecular_mass(self, molecule):
        """Get molecular mass for common molecules"""
        masses = {
            '14N-1H3': 17.0 * u.u,  # NH3
            '1H2-16O': 18.0 * u.u,  # H2O  
            '12C-1H4': 16.0 * u.u,  # CH4
        }
        return masses.get(molecule, 20.0 * u.u)  # Default fallback