# Physics Background

This page describes the physical models and equations implemented in XRayLabTool.

## Refractive index for X-rays

For X-rays, the complex refractive index of a material is:

```math
n = 1 - \delta - i\beta
```

where δ (dispersion) and β (absorption) are both small positive numbers (typically 10⁻⁵ to 10⁻⁸).

## Optical constants δ and β

The dispersion and absorption coefficients are calculated from:

```math
\delta = \frac{\lambda^2}{2\pi} r_e \frac{10^{6}\,\rho N_A}{M} \sum_i n_i f_{1,i}(E)
```

```math
\beta = \frac{\lambda^2}{2\pi} r_e \frac{10^{6}\,\rho N_A}{M} \sum_i n_i f_{2,i}(E)
```

where:
- ``\lambda`` — X-ray wavelength (m)
- ``r_e = 2.818 \times 10^{-15}`` m — classical electron radius (Thomson scattering length)
- ``\rho`` — mass density (g/cm³); the factor ``10^6`` converts it to g/m³ so that ``\lambda`` in metres gives a dimensionless result
- ``N_A`` — Avogadro's number
- ``M`` — molecular weight (g/mol)
- ``n_i`` — number of atoms of element ``i`` per formula unit
- ``f_{1,i}(E), f_{2,i}(E)`` — energy-dependent atomic scattering factors for element ``i``

## Atomic scattering factors

The scattering factors ``f_1`` and ``f_2`` describe how an atom scatters X-rays. They are tabulated as functions of photon energy and show sharp features at absorption edges.

At high energies far from absorption edges, ``f_1 \to Z`` (the atomic number), meaning the atom scatters like ``Z`` free electrons.

The data tables are from the CXRO (Center for X-ray Optics) Henke tables, covering energies from 30 eV to 30 keV for elements H through U.

## Derived quantities

### Energy–wavelength conversion

```math
\lambda [\text{Å}] = \frac{12.39842}{E [\text{keV}]}
```

### Critical angle for total external reflection

```math
\theta_c = \sqrt{2\delta} \quad (\text{radians})
```

### Attenuation length

The 1/e intensity attenuation length:

```math
L = \frac{\lambda}{4\pi\beta}
```

### Scattering length density (SLD)

```math
\text{SLD} = \frac{2\pi}{\lambda^2} (\delta + i\beta)
```

SLD is energy-independent to first order (the ``\lambda^2`` in δ and β cancels with the ``1/\lambda^2`` prefactor), making it useful for reflectometry.

### Electron density

```math
\rho_e = \frac{\rho N_A Z_{mol}}{M} \times 10^{-24}
```

where ``Z_{mol}`` is the total number of electrons per formula unit. With ``\rho`` in g/cm³ and ``M`` in g/mol, the factor ``10^{-24}`` (cm³ → Å³) gives ``\rho_e`` in electrons/Å³.

## Physical constants

The package uses the following constants (as defined in `src/XRayLabTool.jl`):

| Constant | Symbol | Value |
|----------|--------|-------|
| Thomson scattering length | ``r_e`` | 2.8179403205 × 10⁻¹⁵ m |
| Speed of light | ``c`` | 299792458 m/s |
| Planck constant | ``h`` | 6.62607015 × 10⁻³⁴ J·s (exact) |
| Elementary charge | ``e`` | 1.602176634 × 10⁻¹⁹ C (exact) |
| Avogadro's number | ``N_A`` | 6.02214076 × 10²³ mol⁻¹ (exact) |

## References

- B.L. Henke, E.M. Gullikson, and J.C. Davis, *X-ray interactions: photoabsorption, scattering, transmission, and reflection at E=50-30000 eV, Z=1-92*, Atomic Data and Nuclear Data Tables **54**, 181-342 (1993).
- [CXRO X-Ray Interactions with Matter](http://henke.lbl.gov/optical_constants/)
- [NIST X-Ray Form Factor, Attenuation, and Scattering Tables](https://www.nist.gov/pml/x-ray-form-factor-attenuation-and-scattering-tables)
