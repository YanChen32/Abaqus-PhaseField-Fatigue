# -*- coding: utf-8 -*-
# Pure-Python helpers for PF-Fatigue v1.0.
# No Abaqus imports are used here, so the set-up logic can be tested outside
# Abaqus/CAE (Python 2.7 and 3.x).

import hashlib
import math
import os

VERSION = '1.0'

# Bundled user subroutines (identical to the files released with the paper).
SUBROUTINES = {
    'cycle': ('PF_Fatigue_Exact_EnergySplit.f',
              '224267f1c4e5cb272e114c19cbb103b821276539bd63f26546d38dd0acb2ff15'),
    'envelope': ('PF_Fatigue_GC_Envelope_Final.f',
                 '6032aa447c9cedce3ab6c2d36a75948c5c3ac3944135fd68753ae46d0840cd65'),
}

SPLITS = (('None', 0), ('Amor', 1), ('Miehe', 2), ('No-tension', 3))

# Element types used in the examples of the paper.
SUPPORTED_ELEMENTS = ('CPE4T', 'CPS4T', 'CPE3T', 'C3D8T')

# Upper bound of the scaled cycle increment of the envelope mode above which Check Model warns.
# The crack-growth calculations of the paper use dN* = 0.625.
DN_STAR_WARN = 1.0

# Base-material choice that lets PF-Fatigue create and assign a solid section itself.
AUTO_SECTION = '(auto) section for all parts'
AUTO_SECTION_NAME = 'PF-Fatigue-Section'

# Components that can drive the fatigue load.
BC_COMPONENTS = ('u1', 'u2', 'u3', 'ur1', 'ur2', 'ur3')
LOAD_COMPONENTS = ('cf1', 'cf2', 'cf3', 'cm1', 'cm2', 'cm3', 'magnitude')


def nonzero_components(obj, names):
    """[(name, value)] of the numeric, non-zero attributes of obj (UNSET/FREED etc. are skipped)."""
    out = []
    for name in names:
        try:
            value = float(getattr(obj, name))
        except Exception:
            continue
        if abs(value) > 1.0e-30:
            out.append((name, value))
    return out


def guess_base_material(model, materialNames):
    """Material used by an existing section; otherwise the automatic-section choice."""
    try:
        for s in model.sections.keys():
            try:
                m = model.sections[s].material
            except Exception:
                continue
            if m in materialNames:
                return m
    except Exception:
        pass
    return AUTO_SECTION


def guess_driver(model):
    """('BC | name', step) of the first Load/BC with exactly one non-zero component, else (None, None)."""
    try:
        steps = [s for s in model.steps.keys() if s != 'Initial']
        groups = (('BC', model.boundaryConditions, BC_COMPONENTS, 'boundaryConditionStates'),
                  ('LOAD', model.loads, LOAD_COMPONENTS, 'loadStates'))
        for kind, repo, comps, attr in groups:
            for n in sorted(repo.keys()):
                for s in steps:
                    try:
                        states = getattr(model.steps[s], attr)
                        if n not in states.keys():
                            continue
                        nz = nonzero_components(states[n], comps)
                    except Exception:
                        continue
                    if len(nz) == 1:
                        return '%s | %s' % (kind, n), s
    except Exception:
        pass
    return None, None


def split_code(name):
    for label, code in SPLITS:
        if label == name:
            return code
    raise ValueError('Unknown energy split: %s' % name)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 16), b''):
            h.update(chunk)
    return h.hexdigest()


def subroutine_path(mode, plugin_dir):
    name, _ = SUBROUTINES[mode]
    return os.path.join(plugin_dir, 'subroutines', name)


def check_subroutine(mode, plugin_dir):
    """Return (errors, notes) for the bundled subroutine of a mode."""
    path = subroutine_path(mode, plugin_dir)
    name, expected = SUBROUTINES[mode]
    if not os.path.isfile(path):
        return ['Bundled user subroutine not found: %s' % path], []
    actual = sha256(path)
    if expected and actual != expected:
        return ['Bundled %s has been modified (SHA256 %s...). Restore the released file.'
                % (name, actual[:12])], []
    return [], ['Bundled subroutine %s verified (SHA256 %s...).' % (name, actual[:12])]


def _num(value, label, errors, kind=float):
    try:
        return kind(value)
    except Exception:
        errors.append('%s must be numeric.' % label)
        return None


def validate_material(E, nu, Gc, ell):
    errors = []
    E = _num(E, 'E', errors); nu = _num(nu, 'nu', errors)
    Gc = _num(Gc, 'Gc', errors); ell = _num(ell, 'l0', errors)
    if E is not None and E <= 0.0:
        errors.append('E must be > 0.')
    if nu is not None and not (-0.999 < nu < 0.4999):
        errors.append('Poisson ratio must satisfy -0.999 < nu < 0.4999.')
    if Gc is not None and Gc <= 0.0:
        errors.append('Gc must be > 0.')
    if ell is not None and ell <= 0.0:
        errors.append('l0 must be > 0.')
    return errors


def alpha_t0(Gc, ell):
    """Reference fatigue energy Gc/(12 ell) used by the envelope law."""
    return float(Gc) / (12.0 * float(ell))


def validate_cycle_inputs(max_value, load_ratio, cycles, ipc, alpha_t, kres):
    errors = []
    mv = _num(max_value, 'Maximum value', errors)
    r = _num(load_ratio, 'Load ratio R', errors)
    n = _num(cycles, 'Number of cycles', errors, int)
    q = _num(ipc, 'IPC', errors, int)
    a = _num(alpha_t, 'alpha_T', errors)
    k = _num(kres, 'k_res', errors)
    if mv is not None and mv == 0.0:
        errors.append('Maximum value must be non-zero.')
    if r is not None and r >= 1.0:
        errors.append('Load ratio R must be less than 1.0.')
    if n is not None and n < 1:
        errors.append('Number of cycles must be at least 1.')
    if q is not None and (q < 4 or q % 4 != 0):
        errors.append('IPC must be a multiple of 4 (>= 4) so the triangular turning points lie on increments.')
    if a is not None and a <= 0.0:
        errors.append('alpha_T must be > 0.')
    if k is not None and k < 0.0:
        errors.append('k_res must be >= 0.')
    return errors


def validate_envelope_inputs(max_value, load_ratio, cycles, delta_n, eta_t, p, kappa):
    errors = []
    mv = _num(max_value, 'Peak value', errors)
    r = _num(load_ratio, 'Load ratio R', errors)
    n = _num(cycles, 'Number of cycles', errors)
    d = _num(delta_n, 'Cycle increment dN', errors)
    e = _num(eta_t, 'eta_T', errors)
    pp = _num(p, 'p', errors)
    k = _num(kappa, 'kappa', errors)
    if mv is not None and mv == 0.0:
        errors.append('Peak value must be non-zero.')
    if r is not None and r >= 1.0:
        errors.append('Load ratio R must be less than 1.0.')
    if n is not None and n <= 0.0:
        errors.append('Number of cycles must be > 0.')
    if d is not None and d <= 0.0:
        errors.append('Cycle increment dN must be > 0.')
    if n is not None and d is not None and d > 0.0 and n / d > 1.0e7:
        errors.append('More than 1e7 increments requested; increase dN or reduce the number of cycles.')
    # the step uses fixed increments of dN, so the step time must be a whole number of increments
    if n is not None and d is not None and n > 0.0 and d > 0.0 and abs(n / d - round(n / d)) > 1.0e-6 * max(1.0, n / d):
        errors.append('Number of cycles must be an integer multiple of dN (fixed increments of dN).')
    if e is not None and e <= 0.0:
        errors.append('eta_T must be > 0.')
    if pp is not None and pp <= 0.0:
        errors.append('p must be > 0.')
    if k is not None and k < 0.0:
        errors.append('kappa must be >= 0.')
    return errors


def scaled_increment(delta_n, load_ratio, eta_t, p, kappa):
    """Scaled cycle increment dN* = dN beta^(2 kappa p) / eta_T, beta = (1-R)/2 (paper Eq. (29))."""
    s = (1.0 - float(load_ratio)) / 2.0
    return float(delta_n) * s ** (2.0 * float(kappa) * float(p)) / float(eta_t)


def cycle_props(E, nu, Gc, ell, alpha_t, kres, split):
    """PROPS of PF_Fatigue_Exact_EnergySplit.f: E, nu, Gc, l, alpha_T, k_res, ISPLIT."""
    return (float(E), float(nu), float(Gc), float(ell), float(alpha_t), float(kres), float(split_code(split)))


def envelope_props(E, nu, Gc, ell, eta_t, split, p, load_ratio, kappa):
    """PROPS of PF_Fatigue_GC_Envelope_Final.f: E, nu, Gc, l, eta_T, ISPLIT, p, R, kappa."""
    return (float(E), float(nu), float(Gc), float(ell), float(eta_t), float(split_code(split)),
            float(p), float(load_ratio), float(kappa))


def build_triangle_amplitude(cycles, load_ratio, ipc):
    """Triangular history starting from zero.

    A quarter-cycle ramp reaches the maximum at t = 0.25; afterwards the load
    cycles between max and R*max; one cycle has pseudo-time 1; dt = 1/IPC.
    Returns (data, time_period, dt, total_increments).
    """
    n = int(cycles)
    r = float(load_ratio)
    q = int(ipc)
    data = [(0.0, 0.0), (0.25, 1.0)]
    for k in range(n):
        data.append((0.75 + float(k), r))
        data.append((1.25 + float(k), 1.0))
    dt = 1.0 / float(q)
    return tuple(data), float(n) + 0.25, dt, n * q + q // 4


def build_envelope_schedule(cycles, delta_n):
    """Constant peak amplitude over the step; step time = cycle count N; dt = dN."""
    n = float(cycles)
    d = float(delta_n)
    total = int(math.ceil(n / d - 1.0e-9))
    return ((0.0, 1.0), (n, 1.0)), n, d, total


def element_size_stats(elements):
    """elements: iterable of lists of nodal coordinate tuples (2D or 3D).

    Returns (h_min, h_median, n) where h is the shortest edge between
    consecutive nodes of each element (corner nodes only).
    """
    sizes = []
    for coords in elements:
        m = len(coords)
        if m < 2:
            continue
        best = None
        for i in range(m):
            for j in range(i + 1, m):
                d = math.sqrt(sum((a - b) ** 2 for a, b in zip(coords[i], coords[j])))
                if d > 0.0 and (best is None or d < best):
                    best = d
        if best is not None:
            sizes.append(best)
    if not sizes:
        return None, None, 0
    sizes.sort()
    return sizes[0], sizes[len(sizes) // 2], len(sizes)


def mesh_advice(ell, h_min, h_median):
    """Notes/warnings on the length-scale resolution in the refined region."""
    notes, warns = [], []
    if h_min is None:
        warns.append('Element size could not be evaluated.')
        return notes, warns
    ratio = float(ell) / h_min
    notes.append('Smallest element edge h_min = %.4g (l0/h_min = %.2f); median edge %.4g.'
                 % (h_min, ratio, h_median))
    if ratio < 2.0:
        warns.append('l0/h_min = %.2f < 2: the crack path region is coarse for the chosen l0. '
                     'Refine the mesh where the crack is expected (h <= l0/2 is recommended).' % ratio)
    return notes, warns


def parse_load_source(text):
    if not text or '|' not in text:
        raise ValueError('Select an existing Load or BC.')
    kind, name = text.split('|', 1)
    kind = kind.strip().upper()
    name = name.strip()
    if kind not in ('BC', 'LOAD'):
        raise ValueError('Unsupported load-source selector.')
    if not name:
        raise ValueError('Load/BC name is empty.')
    return kind, name
