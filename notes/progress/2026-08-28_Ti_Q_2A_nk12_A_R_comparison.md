# Ti_Q_2A 12x12x12: in-house MMN A(R) versus postw90 transl_inv_full

## Question

Do the in-house and postw90 position-matrix constructions use compatible phase
and real-space conventions, and does either choice explain the curvature
difference between trial 03 (48 WFs) and trial 04 (38 WFs)?

## Finite-difference phases

For a Wannier-gauge overlap `M_ij(k,b)`, the in-house construction forms

```text
A^ours_a,ij(k) = i sum_b w_b b_a M_ij(k,b) exp(i b.tau_j),
```

Fourier transforms the completed `A(k)`, and then distributes each orbital-pair
block among the images minimizing

```text
|R + tau_j - tau_i|.
```

With `transl_inv_full=true`, postw90 instead applies the symmetric link phase

```text
phase1 = exp[i b.(tau_i + tau_j)/2]
```

before Fourier transformation.  Inside the per-link Wigner--Seitz mapping it
then applies

```text
phase2 = exp[-i b.R_ij/2],
```

where `R_ij` is the image selected for that orbital pair.  The two schemes are
therefore different finite-mesh approximations, not differently named copies
of the same `A(R)`.  They should approach one another as the finite-difference
links become shorter.

Both schemes use the same mathematical real-space choice, minimizing
`|R + tau_j - tau_i|` and averaging tied representatives.  In the postw90
route, postw90 maps `A(R)` internally and the repository maps `H(R)` with that
same rule; the exported `A(R)` is not mapped a second time.

## R support correction

Both 12x12x12 constructions contain exactly the same union of 2,197 R vectors
(`[-6,6]^3`), with no vector unique to either construction.  Trial 04's mapped
`H(R)` contains only 1,495 nonzero blocks.  The old MMN installer consequently
dropped 702 `A(R)` blocks that were absent from the nonzero H support, although
their combined relative A norm was only `5.063e-6`.

The installer now keeps the full, independent A support and adds explicit zero
H blocks, matching the existing postw90 installer.  This changes the 12x12x12
trial comparison negligibly: the combined relative curve difference moves from
1.7256% to 1.7262%.

## Direct A(R) comparison

After placing both matrices in the common pair-dependent R convention:

| trial | relative Frobenius difference in A(R) | largest element difference |
|---|---:|---:|
| 03 (48 WFs) | 0.07134% | 0.001088 Ang |
| 04 (38 WFs) | 0.12192% | 0.001477 Ang |

## Same-H curvature comparison: change only A(R)

The repository curvature evaluator was run at the exact same 202 postw90 path
points with the same mapped `H(R)`.  Only the installed position matrix was
changed.

| trial | normalized all-component curve difference | largest absolute change |
|---|---:|---:|
| 03 (48 WFs) | 0.3257% | 0.008193 Ang^2 |
| 04 (38 WFs) | 0.5320% | 0.008574 Ang^2 |

The small `xy` curve has larger component-relative changes (3.65% and 15.5%),
but the absolute changes are only 0.00216 and 0.00772 Ang^2.

## Same-A evaluator validation

Supplying postw90's exported `A(R)` to the repository evaluator gives:

| trial | repository versus postw90 curve difference | largest absolute difference |
|---|---:|---:|
| 03 | 0.01976% | 0.0003560 Ang^2 |
| 04 | 0.00000683% | 4.09e-7 Ang^2 |

Thus the repository evaluator reproduces postw90 when both use the same
position matrix.  The larger, but still small, trial-03 residual is consistent
with the lower precision of the printed `_hr.dat` used by the repository.

## Trial-set dependence on the generic line

| position construction | all-component difference | max absolute difference |
|---|---:|---:|
| postw90 `transl_inv_full` | 1.6176% | 0.02110 Ang^2 |
| in-house MMN | 1.7262% | 0.02205 Ang^2 |

For the in-house construction the component differences are 1.756% (`yz`),
1.423% (`zx`), and 24.87% (`xy`).  For postw90 they are 1.471%, 1.688%, and
12.07%, respectively.  The relative `xy` values are amplified because `xy` is
the smallest component.

## Source nodes and symmetry

At all 1,728 original 12x12x12 nodes, the in-house curvatures from the two
trials differ by 5.765% in the combined norm (maximum 0.2185 Ang^2), even though
their ab-initio occupied projectors agree to `3.32e-14` in Frobenius norm.  The
equality of projectors at isolated nodes does not fix their derivatives or the
finite-difference connection.

Each interpolant separately passes the surviving off-grid symmetries.  Across
24 generic points, the largest absolute residuals are:

| trial | T | C4z | Mx |
|---|---:|---:|---:|
| 03 | 8.30e-10 | 1.05e-8 | 9.62e-9 Ang^2 |
| 04 | 7.75e-10 | 1.24e-9 | 1.08e-9 Ang^2 |

## Conclusion

There is no evidence here for an R-set mismatch or a bug in the repository
curvature evaluator.  The two A constructions are close but genuinely
different at finite mesh density, and both preserve symmetry.  Both also retain
a similar trial-set dependence.  The construction difference decreases from
8x8x8 to 12x12x12, but the postw90 trial-set difference is not monotonically
smaller (about 1.52% at 8x8x8 and 1.62% at 12x12x12).  A third source mesh is
needed before claiming numerical convergence of the trial independence.

Machine-readable results are in:

- `calculations/SrTiO3/Ti_Q_2A/postw90_reference_comparison_nk12/`
- `calculations/SrTiO3/Ti_Q_2A/curvature_trial_comparison_nk12_full_A_support/`
- `calculations/SrTiO3/Ti_Q_2A/trial_curvature_geometry_nk12/`
