# PythTB issue draft: `model.lattice.orb_vecs = ...` is a silent no-op

**Target:** github.com/pythtb/pythtb (the `external` branch API), issue.

---

**Title:** Assigning to `TBModel.lattice.orb_vecs` silently changes a copy

**Summary.** `TBModel.lattice` returns `copy.copy(self._lattice)` (`pythtb/tbmodel.py:233`), so
`model.lattice.orb_vecs = tau` modifies a throwaway object and the model keeps its old orbital positions. Nothing
warns.

**Effect in practice.** A user who re-embeds a Wannier model this way gets H(k) and velocities with the old positions
while code that reads the `tau` variable uses the new ones. In our Berry-curvature work this mixed two embeddings and
gave a spurious term linear in τ; a finite difference between two models then contained an embedding change.

**Suggestion.** Make `orb_vecs` on the returned lattice read-only (or have `lattice` return a read-only view), or
raise/warn on assignment, and point to `W90.model(orb_vecs=...)` / the constructor as the way to set positions.
