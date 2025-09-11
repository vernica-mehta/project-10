import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.constants as c
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import solve_ivp, cumulative_trapezoid
import opac
from scipy.special import expn
import sys, time
import utils

class GreyModel(object):

    def __init__(self, Teff=124.4, Tirr=6500, g=2288, r=c.R_sun, D=7.78*1e12):
        # setting quantities with units
        self.Teff = Teff * u.K
        self.Tirr = Tirr * u.K
        self._g = g * u.cm / u.s**2
        self.r = r
        self.D = D * u.m

        # initial tau grid
        self._tau_h = np.logspace(-6,6,50)

        # frequency grid
        wavs = np.linspace(1,200,1000) * u.um
        self.freqs = (c.c / wavs).to(u.Hz)

        # cache for opacities and spectrum
        self._opacities = None
        self._spectrum = None

    @property
    def g(self):
        return self._g

    @property
    def tau_h(self):
        return self._tau_h

    @property
    def T_tau(self):
        return utils.T_tau(self._tau_h, self.Teff, self.Tirr, 1, 1, 1/2, self.r, self.D)

    def load_opacities(self):
        import sys, time
        print("[GreyModel] load_opacities: Starting...", file=sys.stderr)
        # Cache results to avoid recomputation
        if self._opacities is not None:
            print("[GreyModel] load_opacities: Using cached result.", file=sys.stderr)
            return self._opacities
        print("[GreyModel] load_opacities: Opening Ross_Planck_opac.fits...", file=sys.stderr)
        t0 = time.time()
        f_kappa = fits.open('Ross_Planck_opac.fits')
        kappa = f_kappa['kappa_Ross [cm**2/g]'].data
        print(f"[GreyModel] load_opacities: Opened Ross_Planck_opac.fits in {time.time()-t0:.2f} s", file=sys.stderr)

        print("[GreyModel] load_opacities: Opening rho_Ui_mu_ns_ne.fits...", file=sys.stderr)
        t1 = time.time()
        f_eos = fits.open('rho_Ui_mu_ns_ne.fits')
        rho = f_eos['rho [g/cm**3]'].data
        print(f"[GreyModel] load_opacities: Opened rho_Ui_mu_ns_ne.fits in {time.time()-t1:.2f} s", file=sys.stderr)

        print("[GreyModel] load_opacities: Setting up grids and interpolators...", file=sys.stderr)
        h = f_kappa[0].header
        T_grid = h['CRVAL1'] + np.arange(h['NAXIS1'])*h['CDELT1']
        Ps_log10 = h['CRVAL2'] + np.arange(h['NAXIS2'])*h['CDELT2']
        # Start ODE at a higher initial pressure to avoid NaN at grid edge
        P0 = 10**(Ps_log10[5])
        f_kappa_interp = RegularGridInterpolator((Ps_log10, T_grid), kappa, bounds_error=False, fill_value=None)
        f_rho_interp = RegularGridInterpolator((Ps_log10, T_grid), rho, bounds_error=False, fill_value=None)

        print("[GreyModel] load_opacities: Solving ODE for pressure profile...", file=sys.stderr)
        t2 = time.time()
        sol = solve_ivp(
            utils.dP_dTau,
            [0, 10**7],
            [P0],
            args=(self.tau_h, self.T_tau, self.g, f_kappa_interp),
            t_eval=self.tau_h,
            method='LSODA'
        )
        print(f"[GreyModel] load_opacities: ODE solve finished in {time.time()-t2:.2f} s", file=sys.stderr)

        Ps = sol.y[0] * u.dyn / u.cm**2
        Ts = self.T_tau
        log10Ps = np.log10(Ps.to_value(u.dyn / u.cm**2))

        Ts_float = Ts.to_value(u.K) if hasattr(Ts, 'unit') else Ts
        interp_points = np.column_stack((log10Ps, Ts_float))
        kappa_bars = f_kappa_interp(interp_points)
        rhos = f_rho_interp(interp_points)

        print("[GreyModel] load_opacities: Done.", file=sys.stderr)
        self._opacities = (Ps, Ts, rhos, kappa_bars)
        return self._opacities

    def apply_opacs(self):

        print("[GreyModel] apply_opacs: Loading opacities...", file=sys.stderr)
        t1 = time.time()
        Ps, Ts, rhos, kappa_bars = self.load_opacities()
        print(f"[GreyModel] apply_opacs: Loaded opacities in {time.time()-t1:.2f} s", file=sys.stderr)

        n_freq = len(self.freqs)
        n_tau = len(self.tau_h)
        log10P_arr = np.log10(Ps.to_value(u.dyn / u.cm**2))
        T_arr = Ts

        print(f"[GreyModel] apply_opacs: Calculating kappa_nu_bars for {n_tau} tau points...", file=sys.stderr)
        # Pass frequency in Hz to kappa_cont, not length
        kappa_nu_bars = np.stack([
            opac.kappa_cont(self.freqs.to_value(u.Hz), log10P_arr[j], T_arr[j]) / rhos[j]
            for j in range(n_tau)
        ], axis=1)

        print(f"[GreyModel] spectrum: Calculating tau_nu for {n_freq} frequencies...", file=sys.stderr)
        t4 = time.time()
        tau_nu_mat = np.array([
            cumulative_trapezoid(kappa_nu_bars[i]/kappa_bars, x=self.tau_h, initial=0)
            for i in range(n_freq)
        ])
        print(f"[GreyModel] spectrum: tau_nu done in {time.time()-t4:.2f} s", file=sys.stderr)

        return T_arr, n_tau, tau_nu_mat

    def make_spec(self, which):

        T_arr, n_tau, tau_nu_mat = self.apply_opacs()

        if which == "internal":
            print(f"[GreyModel] spectrum: Calculating Planck function for {n_tau} tau points...", file=sys.stderr)
            spec = np.stack([
                utils.planck(self.freqs, T_arr[j])
                for j in range(n_tau)
            ], axis=1)
        elif which == "incoming":
            print(f"[GreyModel] spectrum: Calculating Irradiation function for {n_tau} tau points...", file=sys.stderr)
            T_irr = np.ones_like(T_arr) * self.Tirr.value
            spec = np.stack([
                utils.irradiation(self.r, self.D, self.freqs, T_irr[j])
                for j in range(n_tau)
            ], axis=1)

        print(f"[GreyModel] spectrum: Computing final vectorized sum...", file=sys.stderr)
        t5 = time.time()
        expn3_0 = expn(3, 0)
        expn4_mat = expn(4, tau_nu_mat)
        dSlambda = spec[:,1:] - spec[:,:-1]
        dtau_nu = tau_nu_mat[:,1:] - tau_nu_mat[:,:-1]
        dexp = expn4_mat[:,:-1] - expn4_mat[:,1:]
        dtau_nu = np.where(dtau_nu == 0, 1e-30, dtau_nu)
        sum_term = np.sum((dSlambda/dtau_nu) * dexp, axis=1)
        self.spectrum = 0.5 * (spec[:,0]*expn3_0 + sum_term)
        print(f"[GreyModel] spectrum: Final sum done in {time.time()-t5:.2f} s", file=sys.stderr)

        return self.spectrum, spec
    
    @property
    def local_spectrum(self):
        return self.make_spec("internal")[1]
    
    @property
    def irradiated_spectrum(self):
        return self.make_spec("incoming")[1]
    
    @property
    def final_spectrum(self):
        return self.local_spectrum + self.irradiated_spectrum
    
if __name__ == "__main__":
    model = GreyModel()
    spec = model.final_spectrum
    local = model.local_spectrum
    irr = model.irradiated_spectrum

    plt.figure(figsize=(10,5))
    plt.plot((c.c/model.freqs).to_value(u.um), spec)
    plt.plot((c.c/model.freqs).to_value(u.um), local, label='Local Emission', linestyle='--')
    plt.plot((c.c/model.freqs).to_value(u.um), irr, label='Irradiation', linestyle=':')
    plt.xlabel('Wavelength (um)')
    plt.ylabel(r'$F_\nu$ (cgs)')
    plt.title('Emergent Spectrum')
    plt.legend()
    plt.grid()
    plt.savefig('spectrum.png', dpi=200, bbox_inches='tight')
    plt.show()