# Biomarkers in the FK torus examples: check against Lorenzo et al. (2024)

Reference: G. Lorenzo et al., *Personalized forecasting of prostate cancer
during active surveillance*, Cancer Res Commun 4(3):617 (2024), Eqs. 5–10.

Code checked:

- body-fitted: `examples/fischerkolmogrov_torus_example_bodyfitted.cpp`, lines 467–497
- immersed: `examples/fischerkolmogrov_torus_example_immersed.cpp`, lines 578–601 (the
  `biomarkers` lambda) and 843–855

## 1. Unknown and tumour region

The equation is solved for **N̂**, the normalized tumour cell density in [0, 1],
whose growth follows the logistic term ρN̂(1−N̂).

Tumour region (paper, p. 623):

    Ω_T = { x ∈ Ω : N̂(x,t) ≥ N̂_th },   N̂_th = 0.15

0.15 is the value of N̂ on the border of the MRI tumour segmentation, so the
same threshold defines the tumour both in the images and in the simulations.

## 2. Definitions and implementation

Both codes use the smoothed indicator of Ω_T

    H(N̂) = 1 / (1 + exp(−k (N̂ − N̂_th))),   k = 200,

so an integral over Ω_T is computed as an integral over Ω weighted by H. In
both codes this is a quadrature sum Σ_q w_q f(N̂(x_q)) H(N̂(x_q)), using the
same quadrature as the weak form: Gauss rules on the body-fitted patch, and
the Algoim quadrature on cut cells plus Gauss rules on interior cells for the
immersed method.

| Paper | Definition | Code (both methods) | csv column | Status |
|---|---|---|---|---|
| Eq. 7 | V_P = ∫_Ω dΩ | ∫_Ω 1 dΩ | `V_P` [mm³] | ✅ |
| Eq. 5 | V_T = ∫_{Ω_T} dΩ | ∫_Ω H dΩ | `int_H` [mm³] | ✅ |
| Eq. 6 | V_N = ∫_{Ω_T} N̂ dΩ | ∫_Ω N̂ H dΩ | (not written) | ✅ |
| Eq. 8 | N̄ = V_N / V_T | `V_N / V_T` | `N_dash` [–] | ✅ |
| Eq. 9 | N_T = N̄ V_T / V_P = V_N / V_P | `N_dash * V_T / V_P` | `N_T` [–] | ✅ |
| Eq. 10 | A_p = ∫_{Ω_T} ρN̂(1−N̂) dΩ | ∫_Ω ρ N̂(1−N̂) H dΩ | `A_p` [mm³/day] | ✅ as written in Eq. 10 (see §4) |
| — | ∫_Ω N̂ dΩ | ∫_Ω N̂ dΩ | `int_N` [mm³] | not a paper biomarker; extra diagnostic |

**Conclusion:** all six biomarkers in both codes are implemented as Eqs. 5–10
are written, and the two codes compute them in the same way.

## 3. Units

- Lengths are computed in box units (1 box unit = L = 50 mm). Volumes are
  multiplied by L³ before writing, so the csv reports mm³.
- Time is not rescaled, so ρ stays in 1/day in both codes. A_p is therefore
  written in mm³/day.
- N̄ and N_T are ratios and need no conversion.
- Check: V_P = 16 655 mm³ in all runs, equal to the exact torus volume
  2π²Rr²L³.

## 4. Open point: how A_p is normalized in the paper's reported values

The code follows Eq. 10 literally: it integrates ρN̂(1−N̂) over the tumour
region Ω_T. The numbers the paper reports for A_p, however, look normalized by
V_T:

- The paper calls A_p the *mean* proliferation activity and reports it in
  **1/day**: the low-risk median is 3.61·10⁻⁴ (range 2.28–4.94·10⁻⁴) 1/day.
  An integral over Ω_T has units of volume/day.
- In the same patient group V_T ranges from 0.01 to 5.09 cc, a factor of about
  500, but A_p varies by only about 2×. An unnormalized integral would scale
  with V_T.
- ρ_median × ⟨N̂(1−N̂)⟩ ≈ 1.86·10⁻³ × 0.2 ≈ 3.7·10⁻⁴ 1/day, which is the
  reported median.

Our runs at t = 100 days (ρ = 5·10⁻³ 1/day):

| | bf r3 | immersed r5 | bf r4 |
|---|---|---|---|
| A_p (csv, Eq. 10) [mm³/day] | 0.218 | 0.215 | 0.209 |
| A_p / V_T [1/day] | 1.045·10⁻³ | 1.054·10⁻³ | 1.066·10⁻³ |

A_p / V_T ≈ ρ × 0.21, on the same scale as the paper's values.

- **For the immersed vs body-fitted comparison this does not matter**: V_T is
  nearly the same in all runs, so both forms of A_p compare the methods equally
  well.
- **It does matter when comparing A_p values with the paper's.** In that case
  use A_p / V_T (and drop the L³ factor). To ask the supervisor: is the paper's
  reported A_p divided by V_T?

## 5. Differences from the paper (not errors)

1. **Smoothed Heaviside.** The paper uses the sharp cut N̂ ≥ 0.15 and does not
   specify how the integral over Ω_T is evaluated. With k = 200, H goes from
   0.12 to 0.88 for 0.14 ≤ N̂ ≤ 0.16. This is our own choice and gives a small
   bias in V_T, V_N and A_p near the tumour edge.
2. **No t = 0 row.** Biomarkers are written only after each time step, from
   t = Δt onward, so the initial condition does not appear in the csv.
3. **Like-for-like comparison.** Both codes evaluate the biomarkers on the
   end-of-step solution, with the same quadrature as the weak form.

## 6. Results of the comparison at t = 100 days (20 steps, Δt = 5 days)

| | bf r3 (h ≈ 1.9 mm) | immersed r5 (h = 1.6 mm) | bf r4 (h ≈ 0.9 mm) |
|---|---|---|---|
| V_P [mm³] | 16 655 | 16 655 | 16 655 |
| V_T [mm³] | 208.7 | 204.1 | 196.3 |
| N̄ [–] | 0.365 | 0.371 | 0.379 |
| N_T [–] | 4.570·10⁻³ | 4.545·10⁻³ | 4.473·10⁻³ |
| A_p [mm³/day] | 0.218 | 0.215 | 0.209 |
| ∫N̂ [mm³] | 89.98 | 89.82 | 89.86 |

The immersed values lie between the two body-fitted resolutions. The
differences come from mesh resolution, not from the method.

Run lines:

```sh
PHYS="-D 1e-3 -p 5e-3 -t 5 -N 20 -R 3 -s 0.75 -a 0.5 -L 50 --xc 7.88 --yc 25 --zc 27.12"
./build/bin/fischerkolmogrov_torus_example_bodyfitted ${=PHYS} -r 3 --refTor 1 --plotPts 200000 -o 00_cmp_bf_r3 --plot
./build/bin/fischerkolmogrov_torus_example_bodyfitted ${=PHYS} -r 4 --refTor 1 --plotPts 200000 -o 00_cmp_bf_r4 --plot
./build/bin/fischerkolmogrov_torus_example_immersed   ${=PHYS} -r 5 -o 00_cmp_im_r5 --plot
```
