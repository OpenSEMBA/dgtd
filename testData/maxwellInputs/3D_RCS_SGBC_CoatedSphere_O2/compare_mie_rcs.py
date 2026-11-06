#!/usr/bin/env python3
"""Overlay the coated-PEC Mie RCS on the coated-sphere DGTD results.

The reference is a PEC core of radius 0.45 m with a 0.05 m coating
(eps_r = 4, mu_r = 2, sigma = 3 S/m). The time factor is exp(-i omega t),
so a passive coating has Im(eps) > 0. Backscatter uses the Bohren-Huffman
normalization

    sigma = (pi / k^2) |sum_n (2n+1) (-1)^n (a_n - b_n)|^2

and a bistatic cut is (4 pi / k^2) |S|^2, with S2 on the E-plane and S1
on the H-plane. Theta = 0 is forward scatter.
"""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import hankel1, hankel2, spherical_jn, spherical_yn

C0 = 299792458.0
EPS0 = 8.8541878128e-12
A_PEC = 0.45
B_COAT = 0.50
EPS_R = 4.0
MU_R = 2.0
SIGMA = 3.0

HERE = Path(__file__).resolve().parent
SUMMARY = (
    HERE.parents[2]
    / "exports/SimulationData/cuda-1/3D_RCS_SGBC_CoatedSphere_O2"
    / "RCSSurface/coated_sphere_ntff/summary"
)


def wiscombe_nstop(size_parameter):
    x = float(np.abs(size_parameter))
    return max(2, int(np.ceil(x + 4.05 * x ** (1.0 / 3.0) + 2.0)))


def spherical_hankel(n, z, kind):
    """Spherical Hankel function h_n^{(kind)} without cancelling j and y."""
    hankel = hankel1 if kind == 1 else hankel2
    return np.sqrt(np.pi / (2.0 * z)) * hankel(n + 0.5, z)


def riccati_outgoing(n, z):
    """xi = z h^{(1)} and its derivative, the wave that decays into a lossy coat."""
    h = spherical_hankel(n, z, 1)
    xi = z * h
    if n == 0:
        xi_p = np.exp(1j * z)
    else:
        h_m = spherical_hankel(n - 1, z, 1)
        xi_p = z * h_m - n * xi / z
    return xi, xi_p


def riccati_incoming(n, z):
    """eta = z h^{(2)} and its derivative."""
    h = spherical_hankel(n, z, 2)
    eta = z * h
    if n == 0:
        eta_p = np.exp(-1j * z)
    else:
        h_m = spherical_hankel(n - 1, z, 2)
        eta_p = z * h_m - n * eta / z
    return eta, eta_p


def riccati(n, z):
    """psi, chi and their derivatives. chi = -z y_n."""
    j = spherical_jn(n, z)
    y = spherical_yn(n, z)
    psi = z * j
    chi = -z * y
    if n == 0:
        psi_p = np.cos(z)
        chi_p = -np.sin(z)
    else:
        j_m = spherical_jn(n - 1, z)
        y_m = spherical_yn(n - 1, z)
        psi_m = z * j_m
        chi_m = -z * y_m
        psi_p = psi_m - n * psi / z
        chi_p = chi_m - n * chi / z
    return psi, chi, psi_p, chi_p


def coating_index(frequency_hz):
    omega = 2.0 * np.pi * frequency_hz
    eps_c = EPS_R + 1j * SIGMA / (omega * EPS0)
    index = np.sqrt(eps_c * MU_R)
    if index.real < 0.0:
        index = -index
    return index


def pec_coated_coefficients(frequency_hz):
    """a_n (TM) and b_n (TE) for a PEC core with one magnetodielectric coat."""
    k = 2.0 * np.pi * frequency_hz / C0
    m = coating_index(frequency_hz)
    k1 = m * k
    x = k * B_COAT
    # Scattered orders are set by the exterior size. Higher orders are
    # negligible and overflow the coating Bessel functions.
    n_stop = wiscombe_nstop(abs(x))
    gamma_te = m / MU_R
    delta_tm = MU_R / m
    electric = np.empty(n_stop, dtype=np.complex128)
    magnetic = np.empty(n_stop, dtype=np.complex128)

    for n in range(1, n_stop + 1):
        xi_a, xi_ap = riccati_outgoing(n, k1 * A_PEC)
        eta_a, eta_ap = riccati_incoming(n, k1 * A_PEC)
        xi_b, xi_bp = riccati_outgoing(n, k1 * B_COAT)
        eta_b, eta_bp = riccati_incoming(n, k1 * B_COAT)
        psi_x, _, psi_xp, _ = riccati(n, x)
        xi, xi_p = riccati_outgoing(n, x)

        # TM: derivative of the coating solution vanishes on the PEC core.
        ratio_tm = -xi_ap / eta_ap
        s_tm = xi_b + ratio_tm * eta_b
        sp_tm = xi_bp + ratio_tm * eta_bp
        electric[n - 1] = (psi_xp * s_tm - delta_tm * sp_tm * psi_x) / (
            xi_p * s_tm - delta_tm * sp_tm * xi
        )

        # TE: the coating solution itself vanishes on the PEC core.
        ratio_te = -xi_a / eta_a
        s_te = xi_b + ratio_te * eta_b
        sp_te = xi_bp + ratio_te * eta_bp
        magnetic[n - 1] = (psi_xp * s_te - gamma_te * sp_te * psi_x) / (
            xi_p * s_te - gamma_te * sp_te * xi
        )
    return electric, magnetic


def bare_pec_coefficients(size_parameter):
    x = float(size_parameter)
    n_stop = wiscombe_nstop(x)
    electric = np.empty(n_stop, dtype=np.complex128)
    magnetic = np.empty(n_stop, dtype=np.complex128)
    for n in range(1, n_stop + 1):
        psi, chi, psi_p, chi_p = riccati(n, x)
        xi = psi - 1j * chi
        xi_p = psi_p - 1j * chi_p
        electric[n - 1] = psi_p / xi_p
        magnetic[n - 1] = psi / xi
    return electric, magnetic


def angular_functions(n_stop, mu):
    """pi_n and tau_n at mu = cos(theta), orders 1..n_stop."""
    pi = np.zeros(n_stop, dtype=np.float64)
    tau = np.zeros(n_stop, dtype=np.float64)
    pi_nm2 = 0.0
    pi_nm1 = 0.0
    for n in range(1, n_stop + 1):
        if n == 1:
            pi_n = 1.0
        else:
            pi_n = ((2.0 * n - 1.0) / (n - 1.0)) * mu * pi_nm1 - (
                n / (n - 1.0)
            ) * pi_nm2
        tau[n - 1] = n * mu * pi_n - (n + 1.0) * pi_nm1
        pi[n - 1] = pi_n
        pi_nm2 = pi_nm1
        pi_nm1 = pi_n
    return pi, tau


def scattering_amplitudes(electric, magnetic, theta):
    """S1 (perpendicular) and S2 (parallel) at scattering angle theta."""
    n_stop = electric.size
    order = np.arange(1, n_stop + 1, dtype=np.float64)
    weight = (2.0 * order + 1.0) / (order * (order + 1.0))
    pi, tau = angular_functions(n_stop, np.cos(theta))
    s1 = np.sum(weight * (electric * pi + magnetic * tau))
    s2 = np.sum(weight * (electric * tau + magnetic * pi))
    return s1, s2


def rcs_from_amplitudes(s1, s2, wavenumber):
    scale = 4.0 * np.pi / wavenumber**2
    return scale * np.abs(s1) ** 2, scale * np.abs(s2) ** 2


def monostatic(electric, magnetic, wavenumber):
    order = np.arange(1, electric.size + 1, dtype=np.float64)
    series = np.sum((2.0 * order + 1.0) * ((-1.0) ** order) * (electric - magnetic))
    return (np.pi / wavenumber**2) * np.abs(series) ** 2


def load_csv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(line for line in handle if not line.lstrip().startswith("#")))


def check_limits():
    """Vacuum coating equals a bare PEC core; a bare PEC matches its own series."""
    frequency = 300e6
    k = 2.0 * np.pi * frequency / C0
    bare = bare_pec_coefficients(k * A_PEC)
    saved = (EPS_R, MU_R, SIGMA)
    globals_ok = True
    # Local override through the function defaults is avoided: call the
    # coefficient routine after temporarily matching a vacuum coating.
    return bare, k, saved, globals_ok


def vacuum_coated_coefficients(frequency_hz):
    """Same solver with eps_r = mu_r = 1 and sigma = 0. Must match PEC(a)."""
    global EPS_R, MU_R, SIGMA
    saved = (EPS_R, MU_R, SIGMA)
    EPS_R, MU_R, SIGMA = 1.0, 1.0, 0.0
    try:
        return pec_coated_coefficients(frequency_hz)
    finally:
        EPS_R, MU_R, SIGMA = saved


def main():
    frequency = 300e6
    k = 2.0 * np.pi * frequency / C0
    coated = vacuum_coated_coefficients(frequency)
    bare = bare_pec_coefficients(k * A_PEC)
    n = min(coated[0].size, bare[0].size)
    err_a = np.max(np.abs(coated[0][:n] - bare[0][:n]))
    err_b = np.max(np.abs(coated[1][:n] - bare[1][:n]))
    if max(err_a, err_b) > 1e-8:
        raise SystemExit(
            f"vacuum coating does not reduce to a PEC core "
            f"(da={err_a:.3e}, db={err_b:.3e})"
        )
    s1, s2 = scattering_amplitudes(*bare, np.pi)
    back_s, _ = rcs_from_amplitudes(s1, s2, k)
    back_series = monostatic(*bare, k)
    if abs(back_s - back_series) / back_series > 1e-8:
        raise SystemExit(
            f"backscatter amplitude {back_s:.6g} disagrees with the series {back_series:.6g}"
        )

    mono = load_csv(SUMMARY / "monostatic_rcs.csv")
    freq = np.array([float(row["frequency_Hz"]) for row in mono])
    dgtd_mono = np.array([float(row["rcs_m2"]) for row in mono])
    mie_mono = np.empty_like(dgtd_mono)
    for i, f_hz in enumerate(freq):
        electric, magnetic = pec_coated_coefficients(f_hz)
        wavenumber = 2.0 * np.pi * f_hz / C0
        mie_mono[i] = monostatic(electric, magnetic, wavenumber)

    fig, ax = plt.subplots()
    ax.plot(freq / 1e6, dgtd_mono, label="DGTD")
    ax.plot(freq / 1e6, mie_mono, "--", label="Mie, coated PEC")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("RCS (m$^2$)")
    ax.set_title("Monostatic RCS, theta = 180 deg")
    ax.legend()
    fig.tight_layout()
    fig.savefig(SUMMARY / "monostatic_rcs_vs_mie.png", dpi=120)
    plt.close(fig)

    for freq_hz, tag in ((200e6, "200MHz"), (400e6, "400MHz")):
        electric, magnetic = pec_coated_coefficients(freq_hz)
        wavenumber = 2.0 * np.pi * freq_hz / C0
        eplane = load_csv(SUMMARY / f"bistatic_{tag}_Eplane.csv")
        hplane = load_csv(SUMMARY / f"bistatic_{tag}_Hplane.csv")
        theta = np.array([float(row["theta_deg"]) for row in eplane]) * np.pi / 180.0
        dgtd_e = np.array([float(row["rcs_m2"]) for row in eplane])
        dgtd_h = np.array([float(row["rcs_m2"]) for row in hplane])
        mie_h = np.empty_like(theta)
        mie_e = np.empty_like(theta)
        for i, th in enumerate(theta):
            s1, s2 = scattering_amplitudes(electric, magnetic, th)
            mie_h[i], mie_e[i] = rcs_from_amplitudes(s1, s2, wavenumber)
        fig, ax = plt.subplots()
        ax.plot(theta * 180.0 / np.pi, dgtd_e, label="DGTD E-plane")
        ax.plot(theta * 180.0 / np.pi, mie_e, "--", label="Mie E-plane")
        ax.plot(theta * 180.0 / np.pi, dgtd_h, label="DGTD H-plane")
        ax.plot(theta * 180.0 / np.pi, mie_h, "--", label="Mie H-plane")
        ax.set_xlabel("Theta (deg)")
        ax.set_ylabel("RCS (m$^2$)")
        ax.set_title(f"Bistatic RCS at {tag}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(SUMMARY / f"bistatic_{tag}_vs_mie.png", dpi=120)
        plt.close(fig)
        back = monostatic(electric, magnetic, wavenumber)
        print(
            f"{tag}: Mie backscatter {back:.4g} m^2, "
            f"DGTD E(180) {dgtd_e[-1]:.4g}, H(180) {dgtd_h[-1]:.4g}"
        )

    rel = np.max(np.abs(mie_mono - dgtd_mono) / np.maximum(mie_mono, 1e-12))
    print(f"monostatic max |DGTD-Mie|/Mie = {rel:.3g}")
    print(f"vacuum-coating coefficient error {max(err_a, err_b):.3e}")

    docs_dir = HERE.parents[2] / "docs"
    pal_mono_path = docs_dir / "export-monostatic.csv"
    if pal_mono_path.is_file():
        pal_mono = load_csv(pal_mono_path)
        f_pal = np.array([float(r["f (MHz)"]) for r in pal_mono])
        rcs_pal_mono = np.array([float(r["RCS (m^2)"]) for r in pal_mono])
        rcs_mie_pec05 = np.array([
            monostatic(*bare_pec_coefficients(2.0 * np.pi * f_hz / C0 * B_COAT), 2.0 * np.pi * f_hz / C0)
            for f_hz in freq
        ])
        fig, axes = plt.subplots(1, 3, figsize=(19, 5.8), dpi=160)
        ax = axes[0]
        ax.plot(freq / 1e6, dgtd_mono, color="#1f77b4", lw=2.2, label="DGTD (SGBC, r = 0.50 m)")
        ax.plot(f_pal, rcs_pal_mono, color="#2ca02c", lw=2.2, marker="o", ms=3.5, markevery=2, label="Palace TF/SF (coated)")
        ax.plot(freq / 1e6, mie_mono, color="#ff7f0e", lw=2.0, ls="--", label="Mie Coated PEC (0.45 / 0.50 m)")
        ax.plot(freq / 1e6, rcs_mie_pec05, color="#d62728", lw=1.8, ls=":", label="Mie Bare PEC (r = 0.50 m)")
        ax.set_title("Monostatic RCS (theta = 180 deg)", fontsize=11, fontweight="bold")
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel("RCS (m$^2$)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8.2)

        for idx, (freq_hz, tag) in enumerate(((200e6, "200MHz"), (400e6, "400MHz"))):
            ax = axes[idx + 1]
            pal_bi = load_csv(docs_dir / f"export-bistatic-{tag}.csv")
            th_e_pal = [float(r["theta (deg)"]) for r in pal_bi if abs(float(r["phi (deg)"])) < 1e-6]
            rcs_e_pal = [float(r["RCS (m^2)"]) for r in pal_bi if abs(float(r["phi (deg)"])) < 1e-6]
            th_h_pal = [float(r["theta (deg)"]) for r in pal_bi if abs(float(r["phi (deg)"]) - 90.0) < 1e-6]
            rcs_h_pal = [float(r["RCS (m^2)"]) for r in pal_bi if abs(float(r["phi (deg)"]) - 90.0) < 1e-6]
            if 0.0 not in th_h_pal and 0.0 in th_e_pal:
                th_h_pal.insert(0, 0.0)
                rcs_h_pal.insert(0, rcs_e_pal[th_e_pal.index(0.0)])
            if 180.0 not in th_h_pal and 180.0 in th_e_pal:
                th_h_pal.append(180.0)
                rcs_h_pal.append(rcs_e_pal[th_e_pal.index(180.0)])

            eplane = load_csv(SUMMARY / f"bistatic_{tag}_Eplane.csv")
            hplane = load_csv(SUMMARY / f"bistatic_{tag}_Hplane.csv")
            th_deg = np.array([float(r["theta_deg"]) for r in eplane])
            th_rad = th_deg * np.pi / 180.0
            dgtd_e = np.array([float(r["rcs_m2"]) for r in eplane])
            dgtd_h = np.array([float(r["rcs_m2"]) for r in hplane])

            k = 2.0 * np.pi * freq_hz / C0
            el_c, mg_c = pec_coated_coefficients(freq_hz)
            el_p, mg_p = bare_pec_coefficients(k * B_COAT)
            mie_c_e, mie_c_h = np.empty_like(th_rad), np.empty_like(th_rad)
            mie_p_e, mie_p_h = np.empty_like(th_rad), np.empty_like(th_rad)
            for i, th in enumerate(th_rad):
                s1, s2 = scattering_amplitudes(el_c, mg_c, th)
                mie_c_h[i], mie_c_e[i] = rcs_from_amplitudes(s1, s2, k)
                s1p, s2p = scattering_amplitudes(el_p, mg_p, th)
                mie_p_h[i], mie_p_e[i] = rcs_from_amplitudes(s1p, s2p, k)

            ax.plot(th_deg, dgtd_e, color="#1f77b4", lw=2.0, ls="-", label="DGTD E-plane")
            ax.plot(th_deg, dgtd_h, color="#1f77b4", lw=2.0, ls="--", label="DGTD H-plane")
            ax.plot(th_e_pal, rcs_e_pal, color="#2ca02c", lw=2.0, ls="-", marker="o", ms=3, markevery=15, label="Palace E-plane")
            ax.plot(th_h_pal, rcs_h_pal, color="#2ca02c", lw=2.0, ls="--", marker="s", ms=3, markevery=15, label="Palace H-plane")
            ax.plot(th_deg, mie_c_e, color="#ff7f0e", lw=1.8, ls=":", label="Mie Coated E-plane")
            ax.plot(th_deg, mie_c_h, color="#ff7f0e", lw=1.8, ls="-.", label="Mie Coated H-plane")
            ax.plot(th_deg, mie_p_e, color="#d62728", lw=1.4, ls="-", alpha=0.65, label="Mie Bare PEC E-plane")
            ax.plot(th_deg, mie_p_h, color="#d62728", lw=1.4, ls="--", alpha=0.65, label="Mie Bare PEC H-plane")
            ax.set_title(f"Bistatic RCS at {tag}", fontsize=11, fontweight="bold")
            ax.set_xlabel("Theta (deg)")
            ax.set_ylabel("RCS (m$^2$)")
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=7.5, ncol=2)

        fig.suptitle("3D Coated Sphere: DGTD vs. Palace vs. Mie (Coated & Bare PEC r=0.50 m)",
                     fontsize=13, fontweight="bold", y=1.01)
        fig.tight_layout()
        fig.savefig(SUMMARY / "comparison_dgtd_palace_mie.png", dpi=160, bbox_inches="tight")
        plt.close(fig)

    print(f"wrote {SUMMARY}")


if __name__ == "__main__":
    main()
