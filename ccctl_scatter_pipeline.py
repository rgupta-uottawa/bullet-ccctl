#!/usr/bin/env python3
"""
ccctl_scatter_pipeline.py
=========================

Full gravitational-lensing pipeline for the Bullet Cluster (1E 0657-56,
z_d = 0.296) in the CCC+TL cosmology, with its direct ΛCDM counterpart.

Computes:
  - Σ_crit for both cosmologies at a chosen reference source redshift
  - Σ_tot(θ) = Σ_gal + Σ_cluster + Σ_gas for both models
       ΛCDM:   Hernquist galaxies + NFW haloes + β-model gas
       CCC+TL: galaxy-scale CCC + cluster-scale CCC + β-model gas (same inputs)
  - κ(θ) = Σ_tot(θ) / Σ_crit(z_s,ref)
  - Deflection α(θ) via DFT convolution
  - Per-image back-projection to the source plane for the 135 gold multiple
    images in bullet_gold_v4.dat, under the convention
         α_phys(z_s) = α_dft × β(z_s) / β(z_s,ref)
    where β(z_s) ≡ D_ds(z_s)/D_s(z_s) is cosmology-specific.
  - System-wise RMS, mean and median source-plane scatter for both models
  - Optional mass-uncertainty band by rerunning with M200 scaled by (1 ± 0.30)

Written in response to Referee Point 1 (deflection normalisation and
source-redshift scaling): the convention above is physically unambiguous
when κ is anchored at a single reference source redshift; the per-source
factor β(z_s)/β(z_s,ref) is identically Σ_crit(z_s,ref)/Σ_crit(z_s).

Usage
-----
  cd code/
  python3 ccctl_scatter_pipeline.py                 # central M200
  python3 ccctl_scatter_pipeline.py --mass-band     # also run ±30% M200

Outputs saved to ../results/.

Requires: numpy, scipy, pandas, matplotlib. Data files in ../data/.

Reproducibility
---------------
Grid: 500×500 over 10 arcmin ⇒ 1.20 arcsec/pixel. Convergence tested against
250×250 (2.40 "/px); aggregate scatter stable to <1.5 % on all metrics.
See docs/CONVENTIONS.md for the exact deflection equation and sign conventions.

Reference parameters
--------------------
    rho_ref × eta_gal   = 4.10×10^-24 g/cm^3   (galaxy turn-off density)
    rho_ref × eta_clus  = 5.61×10^-26 g/cm^3   (cluster-scale turn-off)
    f_star  = 0.015,  f_gas = 0.12   (Paraficz et al. 2016)
    M200_main = 1.5×10^15 M_sun,  M200_bul = 1.5×10^14 M_sun
                                   (Springel & Farrar 2007; Randall et al. 2008)
    NFW c:   main = 1.94,  bullet = 7.12
    Gas:  dual β-model, β=0.65; main (S0=3e8, rc=278 kpc, q=0.75) +
          bullet (S0=8e8, rc=65 kpc, q=0.85) following Brownstein & Moffat 2007

Author: R. P. Gupta & N. Samaras
Last updated: 2026-10-04
"""
import argparse, json, os, time, warnings
import numpy as np, pandas as pd
from scipy.integrate import quad
from scipy.interpolate import CubicSpline, interp1d, RegularGridInterpolator
from scipy.special import gamma as spgamma
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT  = os.path.join(HERE, '..', 'results')
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- Constants
kpc_cm, Msun_g           = 3.0857e21, 1.989e33
Mpc_m, kpc_m, Msun       = 3.0857e22, 3.0857e19, 1.989e30
G_SI, c_SI               = 6.674e-11, 2.998e8
Omega_m, Omega_L, H0_km  = 0.30, 0.70, 70.0
Z_CLUSTER, Z_SOURCE_REF  = 0.296, 2.0

rho_ref_cgs = 2.0e-24                     # g/cm^3 (times ETA_* to give physical product)
rho_ref     = rho_ref_cgs*kpc_cm**3/Msun_g
ETA_GAL     = 2.05                        # eta_gal * rho_ref = 4.10e-24 g/cm^3
ETA_CLUS    = 0.028050                    # eta_clus * rho_ref = 5.61e-26 g/cm^3
# Equivalent to paper's (eta_gal=0.976, eta_clus=0.0134) at rho_ref=4.2e-24.

RA_CEN, DEC_CEN = 104.6587173, -55.9571513
GRID_ARCMIN, NPIX = 10.0, 500

# ---------------------------------------------------------------- Cosmology
def E(z):
    return np.sqrt(Omega_m*(1+z)**3 + Omega_L)

def dc(z):
    v, _ = quad(lambda zp: 1/E(zp), 0, z, limit=200)
    return (c_SI*1e-3/H0_km)*v

def DA_LC(z1, z2):
    return (dc(z2) - dc(z1))/(1 + z2)

_ctl = pd.read_csv(os.path.join(DATA, 'CTL_zobs_DA.csv'))
_spCT = CubicSpline(_ctl['z_obs'].values, _ctl['DA Gpc'].values*1e3)

def make_cosm(label, z_s=Z_SOURCE_REF):
    """Return (kpc_per_arcsec_at_lens, Sigma_crit_at_z_s, Dd, Ds, Dds)."""
    if label == 'LCDM':
        Dd = DA_LC(0, Z_CLUSTER); Ds = DA_LC(0, z_s); Dds = DA_LC(Z_CLUSTER, z_s)
    else:
        Dd = float(_spCT(Z_CLUSTER)); Ds = float(_spCT(z_s))
        Dds = ((1 + z_s)*Ds - (1 + Z_CLUSTER)*Dd)/(1 + z_s)
    kpc_as = Dd*1e3/206265.0
    Sc = (c_SI**2/(4*np.pi*G_SI))*(Ds*Mpc_m/(Dd*Mpc_m*Dds*Mpc_m))*kpc_m**2/Msun
    return kpc_as, Sc, Dd, Ds, Dds

def beta_cosm(z_s, label):
    """D_ds(z_s)/D_s(z_s) in the given cosmology."""
    if label == 'LCDM':
        Ds = DA_LC(0, z_s); Dds = DA_LC(Z_CLUSTER, z_s)
    else:
        Ds = float(_spCT(z_s))
        Dds = ((1 + z_s)*Ds - (1 + Z_CLUSTER)*float(_spCT(Z_CLUSTER)))/(1 + z_s)
    return Dds/Ds

# ---------------------------------------------------------------- Profiles
BETA_GAS = 0.65

def M_gas_3D(S0, rc, r_max):
    rho0 = S0*spgamma(3*BETA_GAS/2)/(np.sqrt(np.pi)*rc*spgamma(3*BETA_GAS/2 - 0.5))
    M, _ = quad(lambda r: 4*np.pi*r**2*rho0*(1 + (r/rc)**2)**(-3*BETA_GAS/2),
                0, r_max, limit=300)
    return M

def rho_eff_gal(r, Ms, ra):
    rho_t = Ms/(16*np.pi*ra**3)
    if r <= ra:
        return (Ms/(2*np.pi))*ra/(r*(r + ra)**3)        # inner Hernquist branch
    u = r/ra
    return 2*np.sqrt(2)*rho_t*(u*(1 + u)**3)**(-0.5)    # outer CCC branch

def rho_eff_clus(r, Ms, ra):
    rho_t = Ms/(8*np.pi*ra**3); u = r/ra
    return 2*np.sqrt(2)*rho_t*(u*(1 + u)**3)**(-0.5)

def sigma_ccc(rho_fn, Ms, ra, n=30):
    R_arr = np.logspace(np.log10(0.01*ra), np.log10(30*ra), n)
    S_arr = np.array([quad(lambda r: 2*rho_fn(r, Ms, ra)*r/np.sqrt(r**2 - R**2),
                           R*(1 + 1e-8), 200*ra, limit=100, epsrel=1e-3)[0]
                      for R in R_arr])
    return interp1d(R_arr, S_arr, bounds_error=False, fill_value=0.)

def sigma_hernquist(R, Ms, ra):
    R = np.atleast_1d(np.float64(R)); out = np.zeros_like(R)
    for i, Ri in enumerate(R):
        s = max(Ri, 1e-5*ra)/ra
        if abs(s - 1) < 1e-4:
            out[i] = (Ms/(2*np.pi*ra**2))/3.
        elif s < 1:
            out[i] = (Ms/(2*np.pi*ra**2*(s**2 - 1)))*(-1 + np.arctanh(np.sqrt(1 - s**2))/np.sqrt(1 - s**2))
        else:
            out[i] = (Ms/(2*np.pi*ra**2*(s**2 - 1)))*(-1 + np.arctan(np.sqrt(s**2 - 1))/np.sqrt(s**2 - 1))
    return np.abs(out)

def nfw_sig(x2d, y2d, xc, yc, r200, cnfw):
    rho_cz = (3*(H0_km*1e3/Mpc_m)**2*E(Z_CLUSTER)**2/(8*np.pi*G_SI))*kpc_m**3/Msun
    r_s = r200/cnfw
    rho_s = (200*rho_cz*cnfw**3)/(3*(np.log(1 + cnfw) - cnfw/(1 + cnfw)))
    R2d = np.sqrt((x2d - xc)**2 + (y2d - yc)**2)
    xa = (R2d/r_s).ravel(); g = np.empty_like(xa)
    m = xa < 1
    if m.any():
        u = xa[m]; g[m] = (1/(u**2 - 1))*(1 - 2/np.sqrt(1 - u**2)*np.arctanh(np.sqrt((1 - u)/(1 + u))))
    m = np.abs(xa - 1) < 1e-4; g[m] = 1/3.
    m = xa > 1
    if m.any():
        u = xa[m]; g[m] = (1/(u**2 - 1))*(1 - 2/np.sqrt(u**2 - 1)*np.arctan(np.sqrt((u - 1)/(u + 1))))
    Sig = (2*rho_s*r_s*g).reshape(x2d.shape); Sig[R2d > r200] = 0.
    return np.abs(Sig)

def sbeta(x2d, y2d, xc, yc, s0, rc, q, pa):
    dx =  (x2d - xc)*np.cos(pa) + (y2d - yc)*np.sin(pa)
    dy = -(x2d - xc)*np.sin(pa) + (y2d - yc)*np.cos(pa)
    return s0*(1 + ((dx**2 + (dy/q)**2)/rc**2))**(-3*BETA_GAS/2 + 0.5)

# ---------------------------------------------------------------- Deflection
def alpha_dft(kap, pix_as):
    """
    Deflection field from the convergence map via DFT convolution.

      alpha(theta) = (1/pi) * integral  kappa(theta') * (theta-theta')/|theta-theta'|^2  d^2 theta'

    kappa is defined at z_s,ref (the Sigma_crit used to build it).  The returned
    alpha_dft(theta) is in arcsec and corresponds to the reduced deflection
    FOR A SOURCE AT z_s,ref.  See 'back_project' for the per-source rescaling.
    """
    N = kap.shape[0]; P = 2*N
    kap_p = np.zeros((P, P)); kap_p[:N, :N] = kap          # zero-pad to avoid wrap-around
    k = np.fft.fftfreq(P)*P
    kx, ky = np.meshgrid(k, k)
    r2 = kx**2 + ky**2; r2[0, 0] = 1.0
    Kx = kx/r2; Ky = ky/r2
    Kx[0, 0] = 0.0; Ky[0, 0] = 0.0
    Kh = np.fft.fft2(kap_p); KxH = np.fft.fft2(Kx); KyH = np.fft.fft2(Ky)
    ax = np.real(np.fft.ifft2(Kh*KxH))[:N, :N]*(pix_as/np.pi)
    ay = np.real(np.fft.ifft2(Kh*KyH))[:N, :N]*(pix_as/np.pi)
    return ax, ay

# ---------------------------------------------------------------- Build Sigma
def build_sigma(fac_main=1.0, fac_bul=None, cache=None, verbose=True):
    """Return (Sigma_LCDM_total, Sigma_CTL_total, cosmology info)."""
    if fac_bul is None: fac_bul = fac_main
    t0 = time.time()

    # Cosmology — reference redshift
    kpc_as_LC,  Sc_LC_ref,  Dd_LC,  Ds_LC_ref,  Dds_LC_ref  = make_cosm('LCDM')
    kpc_as_CTL, Sc_CTL_ref, Dd_CTL, Ds_CTL_ref, Dds_CTL_ref = make_cosm('CTL')
    f_scale = kpc_as_CTL/kpc_as_LC

    # Catalogue
    raw = pd.read_csv(os.path.join(DATA, 'cluster_members_specz_cat.dat'),
                      sep=r'\s+', comment='#',
                      names=['id', 'ra', 'dec', 'mag_f277w', 'specz', 'specz_ref'])
    dra_as  = (raw['ra']  - RA_CEN)*3600*np.cos(np.radians(DEC_CEN))
    ddec_as = (raw['dec'] - DEC_CEN)*3600
    raw['Mstar'] = 10**(-0.406*(raw['mag_f277w'] - 20) + 9.954)
    raw['x_LC']  = dra_as*kpc_as_LC;   raw['y_LC']  = ddec_as*kpc_as_LC
    raw['x_CTL'] = dra_as*kpc_as_CTL;  raw['y_CTL'] = ddec_as*kpc_as_CTL

    # Cluster-scale physical quantities
    M200m = 1.5e15*fac_main; M200b = 1.5e14*fac_bul
    r200m = 2136.0*fac_main**(1/3); r200b = 995.0*fac_bul**(1/3)
    Mgasm, Mgasb   = 0.12*M200m, 0.12*M200b
    Mstarm, Mstarb = 0.015*M200m, 0.015*M200b
    S0m = 3.0e8*(Mgasm/M_gas_3D(3.0e8, 278.0, r200m))
    S0b = 8.0e8*(Mgasb/M_gas_3D(8.0e8,  65.0, r200b))

    pix_as = GRID_ARCMIN*60/NPIX
    lin_as = np.linspace(-GRID_ARCMIN*60/2, GRID_ARCMIN*60/2, NPIX)
    th_x, th_y = np.meshgrid(lin_as, lin_as)
    x2_LC,  y2_LC  = th_x*kpc_as_LC,  th_y*kpc_as_LC
    x2_CTL, y2_CTL = th_x*kpc_as_CTL, th_y*kpc_as_CTL

    X_bul_LC, Y_bul_LC   = -670.0, 193.0
    X_bul_CTL, Y_bul_CTL = X_bul_LC*f_scale, Y_bul_LC*f_scale

    # Gas
    Sgas_LC  = (sbeta(x2_LC,  y2_LC,  -100.,       0., S0m, 278.,  0.75, np.radians(90)) +
                sbeta(x2_LC,  y2_LC,  -570.,     150., S0b,  65.,  0.85, np.radians(15)))
    Sgas_CTL = (sbeta(x2_CTL, y2_CTL, -100.*f_scale,   0.,         S0m, 278.*f_scale, 0.75, np.radians(90)) +
                sbeta(x2_CTL, y2_CTL, -570.*f_scale, 150.*f_scale, S0b,  65.*f_scale, 0.85, np.radians(15)))

    # NFW DM (LCDM only)
    Sdm_LC = (nfw_sig(x2_LC, y2_LC,       0.,       0., r200m, 1.94) +
              nfw_sig(x2_LC, y2_LC, X_bul_LC, Y_bul_LC, r200b, 7.12))

    # Cluster-scale CCC (CTL only)
    ra_cm = (Mstarm/(8*np.pi*ETA_CLUS*rho_ref))**(1/3)
    ra_cb = (Mstarb/(8*np.pi*ETA_CLUS*rho_ref))**(1/3)
    R_main_CTL = np.sqrt(x2_CTL**2 + y2_CTL**2)
    R_bul_CTL  = np.sqrt((x2_CTL - X_bul_CTL)**2 + (y2_CTL - Y_bul_CTL)**2)
    Sclus_CTL = (sigma_ccc(rho_eff_clus, Mstarm, ra_cm)(R_main_CTL) +
                 sigma_ccc(rho_eff_clus, Mstarb, ra_cb)(R_bul_CTL))

    # Galaxy maps — these depend ONLY on individual M* from catalogue, not M200.
    # Use cache if available to avoid rebuilding on mass-band sweeps.
    if cache is None or 'Sbar_LC' not in cache:
        if verbose: print(f"  Building galaxy maps (219 galaxies)...")
        Sbar_LC  = np.zeros((NPIX, NPIX))
        Sgal_CTL = np.zeros((NPIX, NPIX))
        for i, gal in raw.iterrows():
            ra_CTL = (gal['Mstar']/(16*np.pi*ETA_GAL*rho_ref))**(1/3)
            ra_LC  = ra_CTL/f_scale
            R2_LC  = np.sqrt((x2_LC  - gal['x_LC'] )**2 + (y2_LC  - gal['y_LC'] )**2)
            R2_CTL = np.sqrt((x2_CTL - gal['x_CTL'])**2 + (y2_CTL - gal['y_CTL'])**2)
            Sbar_LC  += sigma_hernquist(R2_LC.ravel(), gal['Mstar'], ra_LC).reshape(NPIX, NPIX)
            Sgal_CTL += sigma_ccc(rho_eff_gal, gal['Mstar'], ra_CTL)(R2_CTL)
        if cache is not None:
            cache['Sbar_LC']  = Sbar_LC
            cache['Sgal_CTL'] = Sgal_CTL
    else:
        Sbar_LC  = cache['Sbar_LC']
        Sgal_CTL = cache['Sgal_CTL']

    S_LC = Sdm_LC + Sgas_LC + Sbar_LC
    S_CT = Sclus_CTL + Sgas_CTL + Sgal_CTL

    cosm = dict(kpc_as_LC=kpc_as_LC, kpc_as_CTL=kpc_as_CTL,
                Sc_LC_ref=Sc_LC_ref, Sc_CTL_ref=Sc_CTL_ref,
                f_scale=f_scale, pix_as=pix_as, lin_as=lin_as, NPIX=NPIX,
                M200m=M200m, M200b=M200b, r200m=r200m, r200b=r200b,
                ra_cm=ra_cm, ra_cb=ra_cb)
    if verbose:
        print(f"  Sigma built  ({time.time() - t0:.1f}s)   "
              f"kappa_peak LC={S_LC.max()/Sc_LC_ref:.3f}  CT={S_CT.max()/Sc_CTL_ref:.3f}")
    return S_LC, S_CT, cosm

# ---------------------------------------------------------------- Back-project
def back_project(S_LC, S_CT, cosm, gold):
    """
    For each image in `gold`, compute source-plane position in both cosmologies
    using convention (a): alpha_phys(z_s) = alpha_dft * beta(z_s)/beta(z_s,ref).

    Returns a DataFrame with columns: sys, id, zs, th_x, th_y,
    bxL, byL (LCDM back-projected source position in arcsec),
    bxC, byC (CTL back-projected).
    """
    kL = S_LC/cosm['Sc_LC_ref']
    kC = S_CT/cosm['Sc_CTL_ref']
    axL, ayL = alpha_dft(kL, cosm['pix_as'])
    axC, ayC = alpha_dft(kC, cosm['pix_as'])
    lin_as = cosm['lin_as']
    iLx = RegularGridInterpolator((lin_as, lin_as), axL.T, bounds_error=False, fill_value=0.)
    iLy = RegularGridInterpolator((lin_as, lin_as), ayL.T, bounds_error=False, fill_value=0.)
    iCx = RegularGridInterpolator((lin_as, lin_as), axC.T, bounds_error=False, fill_value=0.)
    iCy = RegularGridInterpolator((lin_as, lin_as), ayC.T, bounds_error=False, fill_value=0.)
    bLCref = beta_cosm(Z_SOURCE_REF, 'LCDM')
    bCTref = beta_cosm(Z_SOURCE_REF, 'CTL')
    rows = []
    for _, img in gold.iterrows():
        pt = np.array([[img['dra_as'], img['ddec_as']]])
        aLx = float(iLx(pt).item()); aLy = float(iLy(pt).item())
        aCx = float(iCx(pt).item()); aCy = float(iCy(pt).item())
        bLC = beta_cosm(img['zs'], 'LCDM'); bCT = beta_cosm(img['zs'], 'CTL')
        rL = bLC/bLCref; rC = bCT/bCTref
        bxL = img['dra_as']  - aLx*rL
        byL = img['ddec_as'] - aLy*rL
        bxC = img['dra_as']  - aCx*rC
        byC = img['ddec_as'] - aCy*rC
        rows.append((img['sys'], img['id'], img['zs'],
                     img['dra_as'], img['ddec_as'],
                     bxL, byL, bxC, byC))
    return pd.DataFrame(rows, columns=['sys', 'id', 'zs', 'th_x', 'th_y',
                                       'bxL', 'byL', 'bxC', 'byC'])

def system_stats(bp, cols):
    sysnames, vals, Ns = [], [], []
    for s, g in bp.groupby('sys'):
        if len(g) < 2: continue
        bx = g[cols[0]].values; by = g[cols[1]].values
        d2 = (bx - bx.mean())**2 + (by - by.mean())**2
        sysnames.append(s); vals.append(np.sqrt(d2.mean())); Ns.append(len(g))
    vals = np.array(vals); Ns = np.array(Ns)
    rms  = np.sqrt(np.sum(vals**2*Ns)/np.sum(Ns))
    mean = vals.mean()
    med  = float(np.median(vals))
    return dict(RMS=float(rms), Mean=float(mean), Median=med,
                N_sys=int(len(sysnames)), systems=sysnames,
                vals=vals.tolist(), Ns=Ns.tolist())

# ---------------------------------------------------------------- Main
def main():
    p = argparse.ArgumentParser(description="CCC+TL Bullet Cluster scatter pipeline")
    p.add_argument('--mass-band', action='store_true',
                   help="Also compute ±30%% M200 variants for the uncertainty band")
    args = p.parse_args()

    print("=" * 72)
    print("CCC+TL Bullet Cluster — source-plane scatter pipeline")
    print("=" * 72)
    print(f"Reference source redshift z_s,ref = {Z_SOURCE_REF}")

    # Load gold images
    gold = pd.read_csv(os.path.join(DATA, 'bullet_gold_v4.dat'), sep=r'\s+',
                       names=['id', 'ra', 'dec', 'sig_ra', 'sig_dec', 'flag', 'zs', 'n'])
    gold['sys']     = gold['id'].str.split('.').str[0]
    gold['dra_as']  = (gold['ra']  - RA_CEN)*3600*np.cos(np.radians(DEC_CEN))
    gold['ddec_as'] = (gold['dec'] - DEC_CEN)*3600
    print(f"Gold images: {len(gold)} images in {gold.sys.nunique()} systems")

    cache = {}
    scenarios = [(1.00, 'central')]
    if args.mass_band:
        scenarios = [(0.70, 'M-30'), (1.00, 'central'), (1.30, 'M+30')]

    results = {}
    for fac, name in scenarios:
        print(f"\n--- Scenario: {name}  (M200 × {fac}) ---")
        S_LC, S_CT, cosm = build_sigma(fac_main=fac, cache=cache)
        bp = back_project(S_LC, S_CT, cosm, gold)
        sL = system_stats(bp, ('bxL', 'byL'))
        sC = system_stats(bp, ('bxC', 'byC'))
        print(f"  Unopt. ΛCDM:  RMS={sL['RMS']:.4f}  Mean={sL['Mean']:.4f}  Median={sL['Median']:.4f}  ({sL['N_sys']} systems)")
        print(f"  CCC+TL:       RMS={sC['RMS']:.4f}  Mean={sC['Mean']:.4f}  Median={sC['Median']:.4f}")
        results[name] = {'fac': fac, 'LCDM': sL, 'CCCTL': sC}
        if name == 'central':
            bp.to_csv(os.path.join(OUT, 'per_image_backprojection.csv'), index=False)

    with open(os.path.join(OUT, 'scatter_results.json'), 'w') as fp:
        json.dump(results, fp, indent=2)
    print(f"\nSaved: results/scatter_results.json, results/per_image_backprojection.csv")

if __name__ == '__main__':
    main()
