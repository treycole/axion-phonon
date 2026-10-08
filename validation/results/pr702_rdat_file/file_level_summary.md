`_r.dat` vs postw90 `write_aa_r` `_r_full.dat` (`transl_inv_full=T`, `use_ws_distance=F`), matched row-for-row

| system | file | WF | nrpts | components | bit-exact vs 6F12.6(ref) | max abs diff | max resid. after rounding | mean rel. diff, |A|>1e-5 | max rel. diff, |A|>1e-5 | Hermiticity (file) |
|---|---|---|---|---|---|---|---|---|---|---|
| SrTiO3 Ti_Q_2A trial-03 | PR #702 _r.dat, 6F12.6 | 48 | 729 | 10,077,696 | 100.000000% | 5.00e-07 | 4.4e-16 | 1.13e-02 | 4.76e-02 | 0.00e+00 |
| SrTiO3 Ti_Q_2A trial-03 | stock _r.dat | 48 | 729 | 10,077,696 | 74.337706% | 2.87e-02 | 2.9e-02 | 4.19e-01 | 2.91e+01 | 1.12e-02 |
| diamond (testw90_example05) | PR #702 _r.dat, 6F12.6 | 4 | 93 | 8,928 | 100.000000% | 4.83e-07 | 1.7e-18 | 1.46e-03 | 5.71e-03 | 0.00e+00 |
| diamond (testw90_example05) | stock _r.dat | 4 | 93 | 8,928 | 63.104839% | 4.38e-02 | 4.4e-02 | 4.99e+00 | 5.42e+01 | 8.61e-02 |
| testw90_basic2 (disentangled) | PR #702 _r.dat, 6F12.6 | 4 | 89 | 8,544 | 100.000000% | 4.99e-07 | 5.6e-17 | 2.85e-04 | 5.98e-03 | 0.00e+00 |
| testw90_basic2 (disentangled) | stock _r.dat | 4 | 89 | 8,544 | 51.264045% | 4.26e-02 | 4.3e-02 | 5.41e+00 | 2.41e+02 | 1.53e-01 |
| diamond, ES24.16E3 rebuild | PR #702 _r.dat, ES24.16E3 | 4 | 93 | 8,928 | n/a (full precision) | 4.44e-16 | n/a (full precision) | 7.37e-14 | 1.32e-12 | 4.59e-16 |

Reference Hermiticity: 4.9e-15.  Rounded-to-zero fraction of the reference under 6F12.6: SrTiO3 Ti_Q_2A trial-03 66.4%, diamond (testw90_example05) 73.1%, testw90_basic2 (disentangled) 57.7%, diamond, ES24.16E3 rebuild 73.1%.
