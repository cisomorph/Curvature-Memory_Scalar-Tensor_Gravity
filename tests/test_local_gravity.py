"""Solar-System formulas must reproduce textbook limits (ROADMAP 1.1)."""
import numpy as np

import audit_driver as drv


def test_brans_dicke_limit():
    """BD (F̃ = Φ, Z = ω/Φ) mapped to a canonical field gives γ − 1 = −1/(ω + 2)."""
    for omega in (10.0, 500.0, 4e4):
        Phi = 1.0
        # canonical ψ with dψ = √Z dΦ  →  dF̃/dψ = 1/√Z;  local_gravity takes F = F̃/2 and F_i = dF/dψ
        F, F_i = Phi / 2, 0.5 / np.sqrt(omega / Phi)
        g1 = drv.local_gravity(F, [F_i])['gamma_minus_1']
        assert abs(g1 - 1 / (omega + 2)) < 1e-12 * max(1.0, 1 / (omega + 2)) + 1e-15


def test_canonical_cmstg_fails_cassini():
    """SIM107 audit: |γ−1| = 9.476e-4 at Ψ̄ = 2.62, Λ0 = 0.003."""
    L0, psi = 0.003, 2.62
    lg = drv.local_gravity(0.5 + L0 * psi**2, [2 * L0 * psi])
    assert abs(lg['gamma_minus_1'] - 9.476e-4) < 2e-6
    assert not lg['cassini_ok']


def test_attractor_switches_off_fifth_force():
    lg = drv.local_gravity(0.5347, [0.0])
    assert lg['gamma_minus_1'] == 0.0 and lg['cassini_ok']
