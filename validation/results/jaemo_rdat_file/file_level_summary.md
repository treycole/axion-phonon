`_r.dat` vs postw90 `write_aa_r` `_r_full.dat` (`transl_inv_full=T`, `use_ws_distance=F`), matched row-for-row

| system | file | WF | nrpts | components | bit-exact vs 6F12.6(ref) | max abs diff | max resid. after rounding | mean rel. diff, |A|>1e-5 | max rel. diff, |A|>1e-5 | Hermiticity (file) |
|---|---|---|---|---|---|---|---|---|---|---|
| diamond_PR702_wsF | PR #702 _r.dat, 6F12.6 | 4 | 93 | 8,928 | 100.000000% | 4.83e-07 | 1.7e-18 | 1.46e-03 | 5.71e-03 | 0.00e+00 |
| diamond_JaeMo_wsF | PR #702 _r.dat, 6F12.6 | 4 | 93 | 8,928 | 100.000000% | 4.83e-07 | 1.7e-18 | 1.46e-03 | 5.71e-03 | 0.00e+00 |
| basic2_PR702_wsF | PR #702 _r.dat, 6F12.6 | 4 | 89 | 8,544 | 100.000000% | 4.99e-07 | 5.6e-17 | 2.85e-04 | 5.98e-03 | 0.00e+00 |
| basic2_JaeMo_wsF | PR #702 _r.dat, 6F12.6 | 4 | 89 | 8,544 | 100.000000% | 4.99e-07 | 5.6e-17 | 2.85e-04 | 5.98e-03 | 0.00e+00 |

Reference Hermiticity: 2.4e-16.  Rounded-to-zero fraction of the reference under 6F12.6: diamond_PR702_wsF 73.1%, diamond_JaeMo_wsF 73.1%, basic2_PR702_wsF 57.7%, basic2_JaeMo_wsF 57.7%.
