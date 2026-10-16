"""Fatigue parameters p, kappa and eta_T of the cycle-domain law from crack-growth data (Section 4.2.2).

The growth rate at (eta_T, kappa, R) is eta_T^-1 S^(2 kappa p) Phi_p(K_max), S = (1 - R)/2, where Phi_p is the rate of
the reference problem (eta_T = 1, R = -1) under the same peak load. Phi_p is obtained from compact-tension (CT)
reference calculations for a few trial values of p; these calculations are not part of this package.

Procedure (Supporting Information, Section 4):
* Paris interval of every measured curve: the longest run of consecutive points with R^2 >= 0.995 and a root-mean-square
  deviation <= 0.035 in log10(da/dN) of the log-log fit, with at least 8 points and at least 40% of the points;
* experimental exponent m: joint fit log10 v = A + m log10 dK + q log10 S over the two calibration curves (R = 0.1, 0.3);
* reference rate: secant rate between successive output frames at the mean crack length of the two frames, after the
  first 20 l0 of crack extension; K_max from the ASTM E647 CT function; two load levels are joined at the
  geometric mean of their overlapping K_max range;
* step 1: p such that the Paris exponent of Phi_p over the K_max range of the calibration data equals m
  (linear interpolation between the bracketing reference calculations, no extrapolation);
* step 2: mean offsets Lambda_j = mean(ln v - ln Phi_p(K_max)) of the two calibration curves give kappa and eta_T.

Input
  paris_experimental_data.csv          material, R, deltaK_MPa_sqrt_m, da_dN_mm_per_cycle
  <ref_dir>/<tag>_REF_S1_eta1_p<100 p>[_P<load>]_..._crack_history.csv
                                       columns N and x_tip_monotonic_mm (crack-tip x coordinate in mm, measured
                                       from the front edge of the specimen); <load> in N/mm, default from GEOMETRY;
                                       written by abaqus_python/extract_CT_crack_history.py
Usage
  python paris_identify.py windows                       Paris intervals and joint fit of both materials
  python paris_identify.py identify <material> <ref_dir> identification and check predictions -> JSON
"""
import csv
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "paris_experimental_data.csv")
CALIBRATION_R = (0.1, 0.3)
# CT models of the article: width W and load-line position x_L (mm), default peak load per unit thickness (N/mm),
# phase-field length scale l0 (mm), file tag and check stress ratios
GEOMETRY = {"AA6005A-T6": dict(W=50.0, xL=12.5, P=412.5, l0=0.1, tag="6005AT6", check=(0.05, 0.5, 0.6)),
            "Q345qC": dict(W=28.0, xL=7.0, P=1030.0, l0=0.1, tag="Q345qC", check=(0.2, 0.4, 0.5))}
TRANSIENT = 20  # excluded crack extension in units of l0


def exp_points(material, R):
    pts = [(float(r["deltaK_MPa_sqrt_m"]), float(r["da_dN_mm_per_cycle"]))
           for r in csv.DictReader(open(DATA, encoding="utf-8-sig"))
           if r["material"] == material and abs(float(r["R"]) - R) < 1e-9]
    return np.array(sorted(p for p in pts if p[0] > 0 and p[1] > 0))


def fit(x, y):
    m, c = np.polyfit(x, y, 1)
    res = y - (m * x + c)
    r2 = 1 - np.sum(res ** 2) / np.sum((y - y.mean()) ** 2)
    return m, c, r2, math.sqrt(np.mean(res ** 2))


def paris_interval(pts):
    """Index range [i, j) of the Paris interval; ties: smaller RMS deviation, larger R^2, lower start."""
    n = len(pts)
    nmin = max(8, math.ceil(0.4 * n))
    x, y = np.log10(pts[:, 0]), np.log10(pts[:, 1])
    best = None
    for i in range(n):
        for j in range(i + nmin, n + 1):
            m, c, r2, rms = fit(x[i:j], y[i:j])
            if r2 >= 0.995 and rms <= 0.035:
                key = (-(j - i), rms, -r2, i)
                if best is None or key < best[0]:
                    best = (key, i, j)
    if best is None:
        raise ValueError("no Paris interval found")
    return best[1], best[2]


def calibration(material):
    curves, info = [], []
    for R in CALIBRATION_R:
        pts = exp_points(material, R)
        i, j = paris_interval(pts)
        curves.append((R, pts[i:j]))
        info.append(dict(R=R, n_total=len(pts), n=j - i, dK=(float(pts[i, 0]), float(pts[j - 1, 0]))))
    rows, rhs = [], []
    for R, pts in curves:
        for dk, v in pts:
            rows.append([1.0, math.log10(dk), math.log10((1 - R) / 2)])
            rhs.append(math.log10(v))
    sol = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)[0]
    return curves, info, dict(A=sol[0], m=sol[1], q=sol[2], n=len(rhs))


def k_ct(alpha):
    return (2 + alpha) / (1 - alpha) ** 1.5 * (0.886 + 4.64 * alpha - 13.32 * alpha ** 2
                                               + 14.72 * alpha ** 3 - 5.6 * alpha ** 4)


def reference_rate(material, path, P):
    """(K_max, v) of one reference calculation after the initial transient, sorted by K_max."""
    g = GEOMETRY[material]
    N, x = [], []
    for r in csv.DictReader(open(path, encoding="utf-8-sig")):
        if r["x_tip_monotonic_mm"]:
            N.append(float(r["N"]))
            x.append(float(r["x_tip_monotonic_mm"]))
    N, x = np.array(N), np.array(x)
    keep = np.r_[True, np.diff(x) > 0]
    N, x = N[keep], x[keep]
    xm = 0.5 * (x[1:] + x[:-1])
    v = np.diff(x) / np.diff(N)
    steady = (xm - x[0]) >= TRANSIENT * g["l0"]
    K = P / math.sqrt(g["W"]) * k_ct((xm - g["xL"]) / g["W"]) / math.sqrt(1000.0)
    o = np.argsort(K[steady])
    return K[steady][o], v[steady][o]


def log_interp(Kq, K, v):
    out = np.exp(np.interp(np.log(Kq), np.log(K), np.log(v)))
    out[(Kq < K.min()) | (Kq > K.max())] = np.nan
    return out


def join_levels(levels):
    """levels: {P: (K, v)}. Joined (K, v): the lower load below, the higher load above the switching K_max."""
    Ps = sorted(levels)
    if len(Ps) == 1:
        return levels[Ps[0]]
    (Kl, vl), (Kh, vh) = levels[Ps[0]], levels[Ps[1]]
    lo, hi = max(Kl.min(), Kh.min()), min(Kl.max(), Kh.max())
    ks = math.sqrt(lo * hi) if lo < hi else Kl.max()
    return np.r_[Kl[Kl < ks], Kh[Kh >= ks]], np.r_[vl[Kl < ks], vh[Kh >= ks]]


def identify(material, ref_dir):
    g = GEOMETRY[material]
    curves, info, jf = calibration(material)
    lo = min(d["dK"][0] / (1 - d["R"]) for d in info)
    hi = max(d["dK"][1] / (1 - d["R"]) for d in info)
    by_p = {}
    for f in glob.glob(os.path.join(ref_dir, g["tag"] + "_REF_S1_eta1_p*_crack_history.csv")):
        b = os.path.basename(f)
        p = float(b.split("_p")[1].split("_")[0]) / 100
        P = float(b.split("_P")[1].split("_")[0]) if "_P" in b.replace("_REF", "") else g["P"]
        by_p.setdefault(p, {})[P] = reference_rate(material, f, P)
    phis, slopes = {}, {}
    for p in sorted(by_p):
        K, v = join_levels(by_p[p])
        phis[p] = (K, v)
        if K.min() <= lo and K.max() >= hi:
            s = (K >= lo) & (K <= hi)
            slopes[p] = fit(np.log(K[s]), np.log(v[s]))[0]
    out = dict(material=material, paris_intervals=info, joint_fit=jf, Kmax_range=(lo, hi),
               reference_exponents={str(p): m for p, m in slopes.items()})
    ps = sorted(slopes)
    ms = [slopes[p] for p in ps]
    if len(ps) < 2 or not min(ms) <= jf["m"] <= max(ms):
        out["status"] = "experimental exponent not bracketed by the reference calculations; add calculations"
        return out
    k = next(i for i in range(len(ps) - 1) if (ms[i] - jf["m"]) * (ms[i + 1] - jf["m"]) <= 0)
    w = (jf["m"] - ms[k]) / (ms[k + 1] - ms[k])
    p = ps[k] + w * (ps[k + 1] - ps[k])

    def phi_p(Kq):
        a, b = log_interp(Kq, *phis[ps[k]]), log_interp(Kq, *phis[ps[k + 1]])
        return np.exp((1 - w) * np.log(a) + w * np.log(b))

    lam = {R: float(np.nanmean(np.log(pts[:, 1]) - np.log(phi_p(pts[:, 0] / (1 - R))))) for R, pts in curves}
    (R1, l1), (R2, l2) = sorted(lam.items())
    S1, S2 = (1 - R1) / 2, (1 - R2) / 2
    kappa = (l1 - l2) / (2 * p * (math.log(S1) - math.log(S2)))
    ln_eta = 2 * kappa * p * math.log(S1) - l1
    out.update(status="ok", p=p, kappa=kappa, eta_T=math.exp(ln_eta), Lambda=lam,
               walker_gamma=2 * kappa * p / jf["m"], predictions={})
    for R in sorted(set(g["check"]) | set(CALIBRATION_R)):
        pts = exp_points(material, R)
        i, j = paris_interval(pts)
        dk, v = pts[i:j, 0], pts[i:j, 1]
        vp = math.exp(-ln_eta) * ((1 - R) / 2) ** (2 * kappa * p) * phi_p(dk / (1 - R))
        ok = np.isfinite(vp)
        ratio = v[ok] / vp[ok]
        out["predictions"]["R=%g" % R] = dict(
            role="calibration" if R in CALIBRATION_R else "check", n=int(ok.sum()), n_interval=j - i,
            within_30_percent=int(np.sum(np.abs(ratio - 1) <= 0.3)),
            slope_exp=fit(np.log(dk[ok]), np.log(v[ok]))[0], slope_model=fit(np.log(dk[ok]), np.log(vp[ok]))[0])
    return out


def main(argv):
    if argv[:1] == ["windows"]:
        for material in GEOMETRY:
            curves, info, jf = calibration(material)
            for d in info:
                print("%-11s R = %.1f  Paris interval dK %.3f-%.3f MPa m^0.5, %d of %d points"
                      % (material, d["R"], d["dK"][0], d["dK"][1], d["n"], d["n_total"]))
            print("%-11s joint fit m = %.4f, q = %.4f" % (material, jf["m"], jf["q"]))
    elif len(argv) == 3 and argv[0] == "identify":
        res = identify(argv[1], argv[2])
        txt = json.dumps(res, indent=2, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
        with open(os.path.join(HERE, "paris_identify_%s.json" % GEOMETRY[argv[1]]["tag"]), "w") as fh:
            fh.write(txt)
        if res["status"] == "ok":
            print("%s: p = %.3f, kappa = %.3f, eta_T = %.2f" % (argv[1], res["p"], res["kappa"], res["eta_T"]))
            for key, d in res["predictions"].items():
                print("    %-7s %-11s within +-30%%: %d of %d" % (key, d["role"], d["within_30_percent"], d["n"]))
        else:
            print(res["status"])
    else:
        print(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
