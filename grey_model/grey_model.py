# grey_model.py
# calculates the planetary spectrum using a grey atmosphere model
# based on code from various sources referenced within and astr4022 class code

#imports
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.constants as c
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.special import expn
import sys
import utils, opac # local files

class GreyModel(object):

    def __init__(self, Teff=144, Tirr=6500, g=2288, r=c.R_sun.value, D=7.78*1e12, verbose=False):
        """Initialise pseudo-grey atmosphere model. Defaults from Jupiter and the Sun.

        PARAMETERS
        ----------
        Teff : `float`
            Effective temperature of the planet in Kelvin. Default is 144 K.
        Tirr : `float`
            Irradiation temperature from the star in Kelvin. Default is 6500 K.
        g : `float`
            Surface gravity of the planet in cm/s^2. Default is 2288 cm/s^2.
        r : `float`
            Radius of the star in meters. Default is the solar radius.
        D : `float`
            Distance from the star to the planet in meters. Default is 7.78e12 m (5.2 AU).
        verbose : `bool`
            If True, self.vprint detailed logs during computation. Default is False.
        """

        # setting quantities with units
        self.Teff = Teff * u.K
        self.Tirr = Tirr * u.K
        self.g = g * u.cm / u.s**2
        self.r = r * u.m
        self.D = D * u.m

        # initial tau grid
        self.tau_h = np.logspace(-6,6,50)

        # frequency and wavelength grids
        wavs = np.linspace(1,200,1000) * u.um
        self.freqs = (c.c / wavs).to(u.Hz)

        # cache for opacities and spectrum
        self._opacities = None
        self._spectrum = None

        self.vprint = print if verbose else lambda *args, **kwargs: None

    @property
    def T_tau(self):
        """Calculate temperature profile as a function of optical depth."""
        return utils.T_tau(self.tau_h, self.Teff, self.Tirr, 1, 1, 1/2, self.r, self.D)

    def load_opacities(self):
        """Load opacity and EOS data, solve for pressure profile, and interpolate opacities.

         RETURNS
         -------
         Ps : `Quantity`
             Pressure profile in dyn/cm^2.
         Ts : `Quantity`
             Temperature profile in K.
         rhos : `ndarray`
             Density profile in g/cm^3.
         kappa_bars : `ndarray`
             Rosseland mean opacity profile in cm^2/g.
        """
        
        self.vprint("[GreyModel] load_opacities: Starting...")

        # Cache results to avoid recomputation
        if self._opacities is not None:
            self.vprint("[GreyModel] load_opacities: Using cached result.")
            return self._opacities
        
        # otherwise, load data and compute
        self.vprint("[GreyModel] load_opacities: Opening Ross_Planck_opac.fits...")
        f_kappa = fits.open('Ross_Planck_opac.fits')
        kappa = f_kappa['kappa_Ross [cm**2/g]'].data # opacities

        self.vprint("[GreyModel] load_opacities: Opening rho_Ui_mu_ns_ne.fits...")
        f_eos = fits.open('rho_Ui_mu_ns_ne.fits')
        rho = f_eos['rho [g/cm**3]'].data # densities

        # interpolations
        self.vprint("[GreyModel] load_opacities: Setting up grids and interpolators...")
        h = f_kappa[0].header
        T_grid = h['CRVAL1'] + np.arange(h['NAXIS1'])*h['CDELT1']
        Ps_log10 = h['CRVAL2'] + np.arange(h['NAXIS2'])*h['CDELT2']
        P0 = 10**(Ps_log10[5])
        f_kappa_interp = RegularGridInterpolator((Ps_log10, T_grid), kappa, bounds_error=False, fill_value=None)
        f_rho_interp = RegularGridInterpolator((Ps_log10, T_grid), rho, bounds_error=False, fill_value=None)

        # solving ODE for pressure profile
        self.vprint("[GreyModel] load_opacities: Solving ODE for pressure profile...")
        sol = solve_ivp(
            utils.dP_dTau,
            [0, 10**7],
            [P0],
            args=(self.tau_h, self.T_tau, self.g, f_kappa_interp),
            t_eval=self.tau_h,
            method='LSODA'
        )

        Ps = sol.y[0] * u.dyn / u.cm**2
        Ts = self.T_tau
        log10Ps = np.log10(Ps.to_value(u.dyn / u.cm**2))

        Ts_float = Ts.to_value(u.K) if hasattr(Ts, 'unit') else Ts
        interp_points = np.column_stack((log10Ps, Ts_float))
        kappa_bars = f_kappa_interp(interp_points)
        rhos = f_rho_interp(interp_points)

        self.vprint("[GreyModel] load_opacities: Done.")
        self._opacities = (Ps, Ts, rhos, kappa_bars)
        return self._opacities

    def apply_opacs(self):
        """Apply opacities to compute frequency-dependent optical depth profiles.
        
         RETURNS
         -------
         T_arr : `Quantity`
             Temperature profile in K.
         n_tau : `int`
             Number of tau points.
         tau_nu_mat : `ndarray`
             Matrix of frequency-dependent optical depths.
        """

        # load opacities and EOS data
        self.vprint("[GreyModel] apply_opacs: Loading opacities...")
        Ps, Ts, rhos, kappa_bars = self.load_opacities()
        n_freq = len(self.freqs)
        n_tau = len(self.tau_h)
        log10P_arr = np.log10(Ps.to_value(u.dyn / u.cm**2))
        T_arr = Ts

        # opacity matrix, using continuum opacities defined in opac.py
        self.vprint(f"[GreyModel] apply_opacs: Calculating kappa_nu_bars for {n_tau} tau points...")
        kappa_nu_bars = np.stack([
            opac.kappa_cont(self.freqs.to_value(u.Hz), log10P_arr[j], T_arr[j]) / rhos[j]
            for j in range(n_tau)
        ], axis=1)

        # optical depth matrix via integration of opacities over initial tau grid
        self.vprint(f"[GreyModel] spectrum: Calculating tau_nu for {n_freq} frequencies...")
        tau_nu_mat = np.array([
            cumulative_trapezoid(kappa_nu_bars[i]/kappa_bars, x=self.tau_h, initial=0)
            for i in range(n_freq)
        ])
        self.vprint(f"[GreyModel] spectrum: tau_nu done")

        return T_arr, n_tau, tau_nu_mat

    def make_spec(self, which):
        """Calculate the emergent spectrum, either local or irradiated.

         PARAMETERS
         ----------
         which : `str`
             "internal" for local emission, "incoming" for irradiation.

         RETURNS
         -------
         spectrum : `ndarray`
             Emergent spectrum as a function of frequency.
        """

        T_arr, n_tau, tau_nu_mat = self.apply_opacs()

        # spectral grids for local and irradiated components
        if which == "internal":
            self.vprint(f"[GreyModel] spectrum: Calculating Planck function for {n_tau} tau points...")
            spec = np.stack([
                utils.planck(self.freqs, T_arr[j])
                for j in range(n_tau)
            ], axis=1)
        elif which == "incoming":
            self.vprint(f"[GreyModel] spectrum: Calculating Irradiation function for {n_tau} tau points...")
            T_irr = np.ones_like(T_arr) * self.Tirr.value / len(T_arr)
            spec = np.stack([
                utils.irradiation(self.r, self.D, self.freqs, T_irr[j])
                for j in range(n_tau)
            ], axis=1)

        # perform vectorized integration on spectral grids using exponential integrals
        self.vprint(f"[GreyModel] spectrum: Computing final vectorized sum...")
        expn3_0 = expn(3, 0)
        expn4_mat = expn(4, tau_nu_mat)
        expn4_mat = np.where(np.isnan(expn4_mat), 0, expn4_mat)
        dSlambda = spec[:,1:] - spec[:,:-1]
        dtau_nu = tau_nu_mat[:,1:] - tau_nu_mat[:,:-1]
        dexp = expn4_mat[:,:-1] - expn4_mat[:,1:]
        dtau_nu = np.where(dtau_nu == 0, 1e-30, dtau_nu)
        sum_term = np.sum((dSlambda/dtau_nu) * dexp, axis=1)
        self.spectrum = 0.5 * (spec[:,0]*expn3_0 + sum_term)

        return self.spectrum
    
    @property
    def local_spectrum(self):
        """Calculate the local emission spectrum."""
        return self.make_spec("internal")
    
    @property
    def irradiated_spectrum(self):
        """Calculate the irradiated spectrum."""
        return self.make_spec("incoming")
    
    @property
    def final_spectrum(self):
        """Calculate the total emergent spectrum."""
        return self.local_spectrum + self.irradiated_spectrum
    
if __name__ == "__main__":
    """Run the grey atmosphere model and plot the emergent spectrum."""

    model = GreyModel()
    spec = model.final_spectrum
    local = model.local_spectrum
    irr = model.irradiated_spectrum

    plt.figure(figsize=(10,5))
    plt.plot((c.c/model.freqs).to_value(u.um), spec, 'k', linewidth=2, label='Total Spectrum')
    plt.plot((c.c/model.freqs).to_value(u.um), local, 'r', linewidth=1, label='Local Emission', linestyle='--')
    plt.plot((c.c/model.freqs).to_value(u.um), irr, 'b', linewidth=1, label='Irradiation', linestyle=':')
    plt.xlabel('Wavelength (um)')
    plt.ylabel(r'$F_\nu$ (cgs)')
    plt.title('Emergent Spectrum')
    plt.legend()
    plt.grid()
    plt.savefig('spectrum.png', dpi=200, bbox_inches='tight')
    plt.show()