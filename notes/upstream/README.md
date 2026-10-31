# Upstream reports (drafts)

Drafts of the bug reports and merge requests for the external codes, one file per bug. **Nothing here has been
sent** except where a file says so. Each draft is written to be pasted as the issue or MR description; the project
story of each bug stays in [`../bugs/`](../bugs/README.md).

| draft | code | kind | project bug | patch / branch | status |
|---|---|---|---|---|---|
| [qe_new_ns_nc_preimage.md](qe_new_ns_nc_preimage.md) | QE | MR + issue | [#7](../bugs/07_qe_new_ns_nc_symmetrizer.md) | `patches/qe_new_ns_nc_preimage.patch`; GitLab `treycole/q-e` `fix/new_ns_nc-preimage` | bug #7 says an MR was submitted 2026-09; **check** it is open, else send |
| [qe_paw_segni_lsign.md](qe_paw_segni_lsign.md) | QE | MR + issue | [#16](../bugs/16_qe_paw_segni_stale_ux.md) | `patches/qe_paw_segni_lsign.patch`; local `fix/paw-segni-lsign-develop` | draft |
| [qe_irreducible_BZ_unfold_trev.md](qe_irreducible_BZ_unfold_trev.md) | QE | MR + issue | [#19](../bugs/19_qe_sym_ftrans_trev_energy.md) | `patches/qe_irreducible_BZ_unfold_trev.patch`; local `fix/irreducible_BZ-unfold-trev` | draft |
| [qe_noncollinear_gga_spurious_state.md](qe_noncollinear_gga_spurious_state.md) | QE | issue (discussion) | [#18](../bugs/18_qe_noncollinear_gga_spurious_state.md) | none | draft |
| [w90_centre_wrap_num_iter0.md](w90_centre_wrap_num_iter0.md) | Wannier90 | issue | [#20](../bugs/20_w90_centre_wrap_long_cell.md) | none | draft |
| [w90_pr702_use_ws_distance.md](w90_pr702_use_ws_distance.md) | Wannier90 | PR comment | [#5](../bugs/05_pr702_rdat_ignores_use_ws_distance.md) | Jae-Mo Lihm's `write_ndegen_applied` branch | draft |
| [pythtb_lattice_copy_noop.md](pythtb_lattice_copy_noop.md) | PythTB | issue | [#2](../bugs/02_pythtb_embedding_noop.md) | none | draft |

**Not needed:** #12 (PAW `v_rad`/`g_rad` uninitialized) is fixed upstream: QEF issue #863, commit `87ac195cf`
(2026-07-10).

**Branches** (`~/Repos/q-e`; all on the GitLab fork `treycole/q-e`, pushed 2026-10-11; MRs not opened yet):
- one per MR, each a single commit on `upstream/develop` (c5a3cf95b, 2026-10-11): `fix/paw-segni-lsign-develop`,
  `fix/irreducible_BZ-unfold-trev`;
- `production/axion-phonon` (fa9b19c8f): the exact source of the production pw.x `db4b99b4` (cluster
  `~/q-e-kfix`) = `fix/new_ns_nc-preimage` (934f8cdd3, on 7.6 develop) + segni (4b480ac68) + irreducible_BZ.

The patches in `patches/` were made with `git format-patch` from those commits; the three QE patches apply cleanly
to upstream develop c5a3cf95b (checked 2026-10-11).
