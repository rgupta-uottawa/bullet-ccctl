# Bullet Cluster in the CCC+TL cosmology — reproducibility release

Companion code and data for

> Gupta, R. P. & Samaras, N., *Gravitational lensing of the Bullet Cluster in
> the CCC+TL cosmology — a controlled proof-of-concept study*, submitted to
> MNRAS (revised, 2026).

## What this release provides

A single, self-contained Python script that reproduces the source-plane
scatter values of Table 4 and Figure 7 from the input catalogues alone, with
every physical assumption and numerical convention stated explicitly.

```
bullet-ccctl-release/
├── README.md                        this file
├── code/
│   └── ccctl_scatter_pipeline.py    the pipeline (~500 lines)
├── data/
│   ├── cluster_members_specz_cat.dat   219 JWST members (Rihtaršič+ 2026)
│   ├── CTL_zobs_DA.csv                 CCC+TL D_A(z) table (Gupta 2024)
│   └── bullet_gold_v4.dat              135 gold multiple images (Rihtaršič+ 2026)
├── docs/
│   └── CONVENTIONS.md               deflection normalisation, in detail
├── figures/
│   └── fig07_source_scatter.pdf     Figure 7 of the paper
└── results/
    ├── scatter_results.json         Table 4 numbers (JSON)
    └── per_image_backprojection.csv per-image back-projection, central M200
```

## Running

```bash
cd code/
python3 ccctl_scatter_pipeline.py                 # central M200
python3 ccctl_scatter_pipeline.py --mass-band     # also ±30% M200 band
```

Requires `numpy`, `scipy`, `pandas`, `matplotlib`. Runtime ≈ 2 minutes on a
modern laptop for the central run; ≈ 2 min 10 s for the ±30% sweep (galaxy
Abel integrals are cached between mass variants).

## Key results

With the convergence anchored at z_s,ref = 2 and the per-source deflection
rescaling α_phys(z_s) = α_dft × β(z_s)/β(z_s,ref) — see `docs/CONVENTIONS.md` —
the pipeline produces, for the 49-system / 135-image gold sample:

| Model                           | k   | RMS (″)              | Mean (″)             | Median (″)           |
|---------------------------------|-----|----------------------|----------------------|----------------------|
| ΛCDM optimised (Lenstool)*      | 32  | 8.28                 | 7.15                 | 8.51                 |
| ΛCDM unoptimised                | 2   | 6.82 (5.75–8.09)     | 5.93 (5.00–7.02)     | 6.71 (5.78–8.27)     |
| CCC+TL                          | 2   | 5.22 (4.63–5.91)     | 4.44 (3.92–5.05)     | 4.92 (4.55–5.59)     |

Ranges in parentheses are the full extent of the M200 = (1 ± 0.30) × 1.5×10¹⁵ M☉
band. *The optimised ΛCDM value is quoted from Rihtaršič et al. (2026); it
cannot be reproduced here without running Lenstool itself.

The CCC+TL unoptimised scatter is smaller than the unoptimised ΛCDM scatter
in the strong-lensing region. Both are substantially smaller than the
optimised Lenstool value on this particular metric, consistent with the
referee's observation that raw source-plane scatter is not a monotonic
goodness-of-fit statistic.

## Physics parameters (as implemented)

| Parameter                   | Value                                   |
|-----------------------------|-----------------------------------------|
| ρ_ref × η_gal               | 4.10 × 10⁻²⁴ g cm⁻³  (galaxy turn-off) |
| ρ_ref × η_clus              | 5.61 × 10⁻²⁶ g cm⁻³  (cluster turn-off)|
| f_⋆  / f_gas                | 0.015 / 0.12                            |
| M200(main) / M200(bul)      | 1.5×10¹⁵ / 1.5×10¹⁴ M☉                 |
| NFW c (main / bul)          | 1.94 / 7.12                             |
| Gas β / r_c (main / bul)    | 0.65 / 278 / 65 kpc                     |
| Σ_crit(ΛCDM, z=2)           | 2.365×10⁹ M☉/kpc²                      |
| Σ_crit(CCC+TL, z=2)         | 1.707×10⁹ M☉/kpc²                      |

The products ρ_ref × η quoted above are equivalently (ρ_ref = 4.2×10⁻²⁴, η_gal
= 0.976, η_clus = 0.0134), as reported in the manuscript; the pipeline uses
the first form because it makes the physical turn-off density manifest.

## Convention statement for the referee

The deflection normalisation and source-redshift scaling used in this
pipeline are stated explicitly in `docs/CONVENTIONS.md`. The convention is
the referee's option (a):

    α_phys(z_s) = α_dft × β(z_s) / β(z_s,ref) ,

equivalently α_phys(z_s) = α_dft × Σ_crit(z_s,ref) / Σ_crit(z_s).

## Data sources

- JWST member catalogue, multiple-image catalogue:
  Rihtaršič et al. 2026, arXiv:2601.22245
- CCC+TL D_A(z) table: Gupta 2024, Universe 10, 266 (appendix)
- NFW main-halo and bullet-halo parameters: Springel & Farrar 2007; Randall
  et al. 2008
- Dual-β gas model normalisation: Brownstein & Moffat 2007
- f_⋆: Paraficz et al. 2016

## Citation

If you use this code or results, please cite the paper above and this release.
