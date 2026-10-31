"""axion_phonon.tube.integrate_at: the weights tile the zone at every (m, planes, half-width), and at the full
half-width it is integrate()'s entry."""
import numpy as np

from axion_phonon import tube


def test_integrate_at_constant_and_matches_integrate():
    kl = tube.klist(nk=8, m=4, half=2, n3=8)
    one = np.ones(len(kl["k"]))
    for m in (1, 2, 4):
        for n3 in (2, 4, 8):
            for half in (0, 1, 2):
                assert abs(tube.integrate_at(one, kl, m, n3, half) - 1.0) < 1e-12, (m, n3, half)
    rng = np.random.default_rng(0)
    c2 = rng.normal(size=len(kl["k"])) + 1j * rng.normal(size=len(kl["k"]))
    levels, n3s, table = tube.integrate(c2, kl)
    for a, m in enumerate(levels):
        for b, n3 in enumerate(n3s):
            assert abs(tube.integrate_at(c2, kl, m, n3, kl["half"]) - table[a, b]) < 1e-12, (m, n3)
