#!/usr/bin/env python3
"""Exact finite-Hilbert-space reference for the PfQMC driven contour."""

import argparse
import json

import numpy as np
from scipy.linalg import eigh, expm


def annihilation(L, i):
    n = 1 << L
    a = np.zeros((n, n), complex)
    for s in range(n):
        if (s >> i) & 1:
            a[s ^ (1 << i), s] = (-1) ** bin(s & ((1 << i) - 1)).count("1")
    return a


def operators(L, delta=1.0, mu=0.0):
    c = [annihilation(L, i) for i in range(L)]
    cd = [x.conj().T for x in c]
    g = [x + y for x, y in zip(c, cd)] + [-1j * (x - y) for x, y in zip(c, cd)]
    A = np.zeros((2 * L, 2 * L), complex)
    for i in range(L - 1):
        for k in range(2):
            a, b = k * L + i, k * L + i + 1
            z = 1j if i % 2 == 0 else -1j
            A[a, b], A[b, a] = z, -z
        A[i, i + 1] += 1j * delta
        A[i + 1, i] += -1j * delta
        A[L + i, L + i + 1] += -1j * delta
        A[L + i + 1, L + i] += 1j * delta
    for i in range(L):
        A[i, L + i] += -1j * mu
        A[L + i, i] += 1j * mu
    K = np.zeros((1 << L, 1 << L), complex)
    for a in range(2 * L):
        for b in range(a + 1, 2 * L):
            K += 0.5 * A[a, b] * (g[a] @ g[b])
    z = [cd[i] @ c[i] - 0.5 * np.eye(1 << L) for i in range(L)]
    zero = np.zeros_like(K)
    even = sum((z[i] @ z[i + 1] for i in range(0, L - 1, 2)), zero.copy())
    odd = sum((z[i] @ z[i + 1] for i in range(1, L - 1, 2)), zero.copy())
    return K, even, odd, z


def integer_slices(value, dt, name):
    x = value / dt
    n = round(x)
    if n < 0 or abs(x - n) > 1e-10:
        raise SystemExit(f"{name}/dt must be a nonnegative integer")
    return n


def slice_operator(Kh, E, O, V, dt):
    # Temporal ket order: K/2 -> V_even -> V_odd -> K/2.
    return Kh @ expm(-dt * V * O) @ expm(-dt * V * E) @ Kh


def normalize_density(rho):
    scale = float(np.trace(rho).real)
    if not np.isfinite(scale) or scale <= 0:
        raise RuntimeError("density matrix acquired a nonpositive/nonfinite trace")
    return rho / scale


def density_observable(rho, op):
    return float((np.trace(rho @ op) / np.trace(rho)).real)


def center_observables(rho, z):
    L = len(z)
    out = []
    for q in (np.pi, np.pi + 2 * np.pi / L):
        rho_q = sum((np.exp(1j * q * i) * z[i] for i in range(L)), np.zeros_like(z[0]))
        out.append(density_observable(rho, rho_q.conj().T @ rho_q) / L**2)
    spi, sdq = out
    return dict(S_pi=spi, S_pi_dq=sdq, R_cdw=1.0 - sdq / spi)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--L", type=int, default=4)
    p.add_argument("--V0", type=float, default=0.0)
    p.add_argument("--Vf", type=float, required=True)
    p.add_argument("--rate", type=float, default=1.0)
    p.add_argument("--dt", type=float, default=0.1)
    p.add_argument("--theta-init", type=float, default=6.0)
    p.add_argument("--beta-trial", type=float, default=8.0)
    p.add_argument("--delta", type=float, default=1.0)
    p.add_argument("--mu", type=float, default=0.0)
    p.add_argument("--ground-diagnostic", action="store_true")
    a = p.parse_args()
    if a.L > 6:
        raise SystemExit("exact reference is restricted to L<=6")
    if a.dt <= 0 or a.beta_trial <= 0:
        raise SystemExit("dt and beta_trial must be positive")

    n_init = integer_slices(a.theta_init, a.dt, "theta_init")
    dv = a.Vf - a.V0
    if abs(dv) < 1e-12:
        n_ramp = 0
    else:
        if a.rate == 0 or dv / a.rate < 0:
            raise SystemExit("rate must point from V0 to Vf")
        n_ramp = integer_slices(dv / a.rate, a.dt, "ramp_time")

    K, E, O, z = operators(a.L, a.delta, a.mu)
    Kh = expm(-0.5 * a.dt * K)

    # This is exactly the projector trial convention.  V_trial=0 makes every
    # full dt slice exp(-dt*K); retain the possible final partial slice.
    rho_trial = np.eye(1 << a.L, dtype=complex)
    rem = a.beta_trial
    n_trial = 0
    while rem > 1e-12:
        step = min(a.dt, rem)
        rho_trial = expm(-step * K) @ rho_trial
        rem -= step
        n_trial += 1
        rho_trial = normalize_density(rho_trial)

    init_slice = slice_operator(Kh, E, O, a.V0, a.dt)
    P0 = np.eye(1 << a.L, dtype=complex)
    for _ in range(n_init):
        P0 = init_slice @ P0
        P0 /= np.linalg.norm(P0)
    rho0 = normalize_density(P0 @ rho_trial @ P0.conj().T)

    schedule = [a.V0 + a.rate * (l + 0.5) * a.dt for l in range(n_ramp)]
    rho = rho0
    for V in schedule:
        U = slice_operator(Kh, E, O, V, a.dt)
        # The right factor is the strict adjoint, hence the reversed bra product.
        rho = normalize_density(U @ rho @ U.conj().T)

    result = dict(
        method="contour_exact_finite_trial_density",
        L=a.L,
        V0=a.V0,
        Vf=a.Vf,
        drive_rate=a.rate,
        dt=a.dt,
        theta_init=a.theta_init,
        beta_trial=a.beta_trial,
        delta=a.delta,
        mu=a.mu,
        boundary="OBC",
        trial_slices=n_trial,
        initial_slices_per_side=n_init,
        n_ramp=n_ramp,
        ramp_time=n_ramp * a.dt,
        ket_slice_bond_order="even_then_odd",
        bra_slice_bond_order="odd_then_even",
        ket_schedule=schedule,
        bra_schedule=list(reversed(schedule)),
        **center_observables(rho, z),
    )

    if a.ground_diagnostic:
        H0 = K + a.V0 * (E + O)
        vals, vecs = eigh(H0)
        keep = np.where(vals - vals[0] < 1e-10)[0]
        rho_gs = vecs[:, keep] @ vecs[:, keep].conj().T
        for V in schedule:
            U = slice_operator(Kh, E, O, V, a.dt)
            rho_gs = normalize_density(U @ rho_gs @ U.conj().T)
        result["ground_subspace_dimension"] = len(keep)
        result["ground_energy"] = float(vals[0])
        result["ground_subspace_reference"] = center_observables(rho_gs, z)

    print(json.dumps(result, separators=(",", ":")))


if __name__ == "__main__":
    main()
