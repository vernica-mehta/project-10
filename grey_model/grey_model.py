import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.constants as c
import astropy.units as u
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import solve_ivp, cumulative_trapezoid
import opac
from scipy.special import expn
import utils


class GreyModel(object):

    def __init__(self, Teff=0, Tirr=6500, g=2288, r=c.R_sun, D=7.78*1e12):
        # setting quantities with units
        self.Teff = Teff * u.K
        self.Tirr = Tirr * u.K
        self._g = g * u.cm / u.s**2
        self.r = r
        self.D = D * u.m

        # initial tau grid
        self._tau_h = np.concatenate((np.arange(3)/3*1e-3,np.logspace(-3,1.3,30)))

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

    @property
    def F_inc(self):
        return utils.irradiation(self.r, self.D, self.freqs, self.Tirr)

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
            [0, 20],
            [P0],
            args=(self.tau_h, self.T_tau, self.g, f_kappa_interp),
            t_eval=self.tau_h,
            method='LSODA'
        )
        print(f"[GreyModel] load_opacities: ODE solve finished in {time.time()-t2:.2f} s", file=sys.stderr)

        Ps = sol.y[0] * u.dyn / u.cm**2
        Ts = self.T_tau
        log10Ps = np.log10(Ps.to_value(u.dyn / u.cm**2))
        # Interpolators expect shape (N,2) for N points
    # Ensure both columns are unitless floats for interpolation
        Ts_float = Ts.to_value(u.K) if hasattr(Ts, 'unit') else Ts
        interp_points = np.column_stack((log10Ps, Ts_float))
        kappa_bars = f_kappa_interp(interp_points)
        # Debug: Check for NaN/inf in kappa_bars and print input/output
        if np.any(np.isnan(kappa_bars)) or np.any(np.isinf(kappa_bars)):
            raise RuntimeError("NaN or inf detected in kappa_bars after interpolation. Aborting.")
        rhos = f_rho_interp(interp_points)
        if np.any(np.isnan(rhos)) or np.any(np.isinf(rhos)):
            raise RuntimeError("NaN or inf detected in rhos after interpolation. Aborting.")

        print("[GreyModel] load_opacities: Done.", file=sys.stderr)
        self._opacities = (Ps, Ts, rhos, kappa_bars)
        return self._opacities

    def spectrum(self, profile=False):
        import sys, time
        if profile:
            import cProfile, pstats, io
            pr = cProfile.Profile()
            pr.enable()
        t0 = time.time()
        print("[GreyModel] spectrum: Starting computation...", file=sys.stderr)
        # Cache the spectrum result to avoid recomputation
        if self._spectrum is not None:
            print("[GreyModel] spectrum: Using cached result.", file=sys.stderr)
            return self._spectrum
        print("[GreyModel] spectrum: Loading opacities...", file=sys.stderr)
        t1 = time.time()
        Ps, Ts, rhos, kappa_bars = self.load_opacities()
        print(f"[GreyModel] spectrum: Loaded opacities in {time.time()-t1:.2f} s", file=sys.stderr)

        n_freq = len(self.freqs)
        n_tau = len(self.tau_h)
        log10P_arr = np.log10(Ps.to_value(u.dyn / u.cm**2))
        T_arr = Ts
        freq_arr = (c.c / self.freqs).value  # in cm

        print(f"[GreyModel] spectrum: Calculating kappa_nu_bars for {n_tau} tau points...", file=sys.stderr)
        t2 = time.time()
        # Pass frequency in Hz to kappa_cont, not length
        kappa_nu_bars = np.stack([
            opac.kappa_cont(self.freqs.to_value(u.Hz), log10P_arr[j], T_arr[j]) / rhos[j]
            for j in range(n_tau)
        ], axis=1)
        
        print(f"[GreyModel] spectrum: kappa_nu_bars done in {time.time()-t2:.2f} s", file=sys.stderr)

        print(f"[GreyModel] spectrum: Calculating Planck function for {n_tau} tau points...", file=sys.stderr)
        t3 = time.time()
        Slambda_mat = np.stack([
            utils.planck(self.freqs, T_arr[j])
            for j in range(n_tau)
        ], axis=1)
        
        print(f"[GreyModel] spectrum: Planck done in {time.time()-t3:.2f} s", file=sys.stderr)

        print(f"[GreyModel] spectrum: Calculating tau_nu for {n_freq} frequencies...", file=sys.stderr)
        t4 = time.time()
        tau_nu_mat = np.array([
            cumulative_trapezoid(kappa_nu_bars[i]/kappa_bars, x=self.tau_h, initial=0)
            for i in range(n_freq)
        ])
        print(f"[GreyModel] spectrum: tau_nu done in {time.time()-t4:.2f} s", file=sys.stderr)

        print(f"[GreyModel] spectrum: Computing final vectorized sum...", file=sys.stderr)
        t5 = time.time()
        expn3_0 = expn(3, 0)
        expn4_mat = expn(4, tau_nu_mat)
        dSlambda = Slambda_mat[:,1:] - Slambda_mat[:,:-1]
        dtau_nu = tau_nu_mat[:,1:] - tau_nu_mat[:,:-1]
        dexp = expn4_mat[:,:-1] - expn4_mat[:,1:]
        dtau_nu = np.where(dtau_nu == 0, 1e-30, dtau_nu)
        sum_term = np.sum((dSlambda/dtau_nu) * dexp, axis=1)
        local = 0.5 * (Slambda_mat[:,0]*expn3_0 + sum_term)
        print(f"[GreyModel] spectrum: Final sum done in {time.time()-t5:.2f} s", file=sys.stderr)

        print(f"[GreyModel] spectrum: Adding irradiation and caching result.", file=sys.stderr)
        irr = (self.F_inc).cgs
        self._spectrum = local + irr

        if profile:
            pr.disable()
            import os
            s = io.StringIO()
            sortby = 'cumulative'
            ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
            ps.print_stats(20)
            with open('spectrum_profile.txt', 'w') as f:
                f.write('PROFILE: spectrum property\n')
                f.write(s.getvalue())

        print(f"[GreyModel] spectrum: Done! Total time: {time.time()-t0:.2f} s", file=sys.stderr)
        return self._spectrum
    
if __name__ == "__main__":
    model = GreyModel()
    spec = model.spectrum(profile=False)
    plt.plot((c.c/model.freqs).to_value(u.um), spec)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Wavelength (um)')
    plt.ylabel(r'$F_\nu$ (cgs)')
    plt.title('Emergent Spectrum')
    plt.grid()
    plt.savefig('spectrum.png', dpi=200, bbox_inches='tight')
    plt.show()