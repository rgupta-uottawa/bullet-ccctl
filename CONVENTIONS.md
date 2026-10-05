# Lensing conventions used in this pipeline

This note documents the normalisation of the deflection field used by
`ccctl_scatter_pipeline.py`, in response to the referee's Point 1.

## Convergence

The convergence map is built once per cosmology at the reference source
redshift z_s,ref = 2.0:

    κ_model(θ) = Σ_model_total(θ) / Σ_crit_model(z_s,ref)

where Σ_crit_model(z_s) depends on the cosmology through its angular-diameter
distances:

    Σ_crit(z_s) = (c² / 4πG) · D_s(z_s) / (D_d · D_ds(z_s)).

At z_d = 0.296 and z_s,ref = 2.0 the pipeline recovers

    Σ_crit_ΛCDM(2)  = 2.365 × 10⁹  M☉ / kpc²
    Σ_crit_CTL(2)   = 1.707 × 10⁹  M☉ / kpc².

(These values match Table 1 of the manuscript.)

## Deflection field

The reduced deflection follows from the standard identity

    α(θ) = (1/π) ∫ κ(θ′) · (θ − θ′) / |θ − θ′|²  d²θ′ ,

computed on the κ-map grid by zero-padded FFT convolution with the kernel

    K(η) = η / |η|²        (η = θ − θ′, in arcsec)

This is the convolution routine in `alpha_dft(kap, pix_as)`. The resulting
α_dft(θ) is the reduced deflection *for a source at z_s = z_s,ref*,
because κ was built with Σ_crit(z_s,ref).

## Per-source rescaling: convention (a)

For an image at a source redshift z_s ≠ z_s,ref, the reduced deflection must
be rescaled. Writing β(z_s) ≡ D_ds(z_s) / D_s(z_s), one has

    Σ_crit(z_s,ref) / Σ_crit(z_s)  =  β(z_s) / β(z_s,ref) .

The per-source physical deflection is therefore

    α_phys(θ; z_s)  =  α_dft(θ)  ·  β(z_s) / β(z_s,ref) .           (★)

Equation (★) is the convention used throughout this pipeline — the referee's
option (a). It is physically unambiguous: κ is anchored at a single reference
redshift and the per-source factor is dimensionless, cosmology-specific, and
tends to unity at z_s = z_s,ref by construction.

The alternative mentioned by the referee,

    α_phys(θ; z_s) = α_dft(θ) · β(z_s)                           [NOT USED]

would be incorrect for a κ anchored at z_s,ref = 2, because it double-counts
the β(z_s,ref) that is already baked into κ via Σ_crit(z_s,ref).

## Source-plane back-projection

With α_phys from (★), the lens equation reads

    β_src = θ − α_phys(θ; z_s) ,

and the per-system scatter is the RMS distance of the images of a given
system from their centroid in the source plane.

## Cosmology dependence of the per-source factor

For the Bullet Cluster multiple-image sample (z_s = 0.9–6.64), the ratio

    β_CTL(z_s) / β_ΛCDM(z_s)

varies only between 1.128 (at z_s = 6.64) and 1.163 (at z_s = 1.60) — a 3%
range. The per-source rescaling therefore cannot by itself produce
large differential effects between the two cosmologies.

## Numerical convergence

The pipeline runs on a 500 × 500 grid over 10′ (1.20 arcsec/pixel). The
aggregate scatter (RMS, Mean, Median) is stable to <1.5% between 500 × 500
and 250 × 250 (2.40 arcsec/pixel) grids. Changing z_s,ref within the sample
range would shift the per-source rescaling pattern but not the physical
deflections, since (★) rebuilds the correct α_phys for any valid reference.
