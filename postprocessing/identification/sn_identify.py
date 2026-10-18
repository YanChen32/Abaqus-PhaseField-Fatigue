"""Fatigue parameters p, eta_T and kappa of the cycle-domain law from stress-life data (Section 4.2.1).

Under homogeneous uniaxial stress the life of the model is

    N_f = eta_T * S^(-2 kappa p) * F_p(r),   S = (1 - R)/2,   r = sigma_max / sigma_c,
    F_p(r) = (6c)^(-p) [ (1 - phi0)^(2p) + c^(-1/2) * int_phi0^phi_f (1 - 4 phi) phi^(-1/2) (1 - phi)^(2p + 1/2) dphi ],

with c = (27/256) r^2, phi0 the smaller root of phi (1 - phi)^3 = c and phi_f = 0.24 (Supporting Information,
Section 2). The integral is evaluated by 128-point Gauss-Legendre quadrature.

Step 1: p is chosen so that the slope of ln F_p(r) against ln r at the stress levels of the R = 0 tests equals the
        slope of ln N_f against ln r of the tests.
Step 2: the mean offsets Lambda_R = mean(ln N_f - ln F_p(r)) of the curves at R = -1 and R = 0 give
        eta_T = exp(Lambda_-1) and kappa = (Lambda_0 - Lambda_-1) / (2 p ln 2).

Input : sn_experimental_data.csv (alloy, R, sigma_a_MPa, N_f_cycles)
Output: sn_identify_result.json and sn_model_curves.csv (model S-N curves for every R)
Usage : python sn_identify.py
"""
import csv
import json
import os

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "sn_experimental_data.csv")
PHI_F = 0.24
# sigma_c: tensile strength used in the article (Table 2); calibration stress ratios R = -1 and 0
ALLOYS = {"AA2024-T3": 490.0, "AA7050-T7451": 524.0}
CALIBRATION_R = (-1.0, 0.0)
GL_X, GL_W = np.polynomial.legendre.leggauss(128)


def phi0(c):
    return brentq(lambda f: f * (1 - f) ** 3 - c, 0.0, 0.25)


def F(r, p, phi_f=PHI_F):
    """Homogeneous life of the reference problem (eta_T = 1, R = -1)."""
    c = 27.0 / 256.0 * r * r
    f0 = phi0(c)
    # substitution u = sqrt(phi): phi^(-1/2) dphi = 2 du
    a, b = np.sqrt(f0), np.sqrt(phi_f)
    u = 0.5 * (b - a) * GL_X + 0.5 * (b + a)
    w = 0.5 * (b - a) * GL_W
    ph = u * u
    integral = np.sum(w * 2 * (1 - 4 * ph) * (1 - ph) ** (2 * p + 0.5))
    return (6 * c) ** (-p) * ((1 - f0) ** (2 * p) + integral / np.sqrt(c))


def S(R):
    return (1 - R) / 2


def load():
    data = {}
    for row in csv.DictReader(open(DATA, encoding="utf-8-sig")):
        data.setdefault(row["alloy"], {}).setdefault(float(row["R"]), []).append(
            (float(row["sigma_a_MPa"]), float(row["N_f_cycles"])))
    return data


def identify(curves, sigma_c):
    """curves: {R: [(sigma_a, N_f)]} containing R = -1 and R = 0."""
    r0 = np.array([sa / (S(0.0) * sigma_c) for sa, _ in curves[0.0]])
    n0 = np.array([n for _, n in curves[0.0]])
    m_exp = -np.polyfit(np.log(r0), np.log(n0), 1)[0]

    def m_model(p):
        return -np.polyfit(np.log(r0), np.log([F(r, p) for r in r0]), 1)[0]

    p = brentq(lambda q: m_model(q) - m_exp, 0.2, 4.0, xtol=1e-12)
    lam = {R: float(np.mean([np.log(n) - np.log(F(sa / (S(R) * sigma_c), p)) for sa, n in pts]))
           for R, pts in curves.items()}
    eta = np.exp(lam[-1.0])
    kappa = (lam[0.0] - lam[-1.0]) / (2 * p * np.log(2))
    return p, eta, kappa, m_exp, lam


def main():
    data = load()
    report = {}
    with open(os.path.join(HERE, "sn_model_curves.csv"), "w", newline="") as fh:
        out = csv.writer(fh)
        out.writerow(["alloy", "R", "sigma_a_MPa", "N_f_cycles"])
        for alloy, sigma_c in ALLOYS.items():
            p, eta, kappa, m_exp, lam = identify({R: data[alloy][R] for R in CALIBRATION_R}, sigma_c)
            res = dict(sigma_c_MPa=sigma_c, p=p, eta_T=eta, kappa=kappa, slope_R0=m_exp,
                       Lambda={str(R): v for R, v in lam.items()},
                       walker_gamma_low_stress=2 * kappa * p / (2 * p + 1), curves={})
            for R, pts in sorted(data[alloy].items()):
                sa = np.array([s for s, _ in pts])
                n_exp = np.array([n for _, n in pts])
                n_mod = np.array([eta * S(R) ** (-2 * kappa * p) * F(s / (S(R) * sigma_c), p) for s in sa])
                ratio = n_exp / n_mod
                res["curves"]["R=%g" % R] = dict(
                    role="calibration" if R in CALIBRATION_R else "check", n=len(pts),
                    within_factor_2=int(np.sum((ratio >= 0.5) & (ratio <= 2.0))),
                    mean_life_ratio_exp_over_model=float(np.exp(np.mean(np.log(ratio)))))
                for r in np.linspace(0.12, 0.985, 120):
                    out.writerow([alloy, R, r * S(R) * sigma_c, eta * S(R) ** (-2 * kappa * p) * F(r, p)])
            report[alloy] = res
            print("%-13s p = %.4f  eta_T = %.1f  kappa = %.4f" % (alloy, p, eta, kappa))
            for key, c in res["curves"].items():
                print("    %-7s %-11s within a factor of 2: %d of %d" % (key, c["role"], c["within_factor_2"], c["n"]))
    with open(os.path.join(HERE, "sn_identify_result.json"), "w") as fh:
        json.dump(report, fh, indent=2)


if __name__ == "__main__":
    main()
