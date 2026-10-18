# -*- coding: utf-8 -*-
"""
Crack initiation angle of the shear-loaded SENT specimen (article Fig. 6(b)).

For each ODB listed in ODB_CASES the script
1. samples about one frame per cycle and takes the first frame in which the crack, traced with phi >= 0.90,
   is at least 0.06 mm long;
2. extracts the crack centerline near the initial tip from NT11 (or NT);
3. fits a straight line to the first 0.03 mm of the centerline (orthogonal regression);
4. reports the angle to the initial horizontal crack for phi = 0.85, 0.90 and 0.95, and their mean and
   standard deviation.

Usage (ODB files in the working directory, names as in ODB_CASES):
    abaqus python extract_SENT2_initiation_angles.py

Output (folder SENT2_initiation_angles in the working directory): one centerline CSV per case and threshold,
SENT2_initiation_angle_summary.csv and a compact table SENT2_initiation_angle_Origin.csv.
"""

from __future__ import print_function
import os
import sys
import csv
import math

import numpy as np
from odbAccess import openOdb


# ======================================================================
# User settings
# ======================================================================

ODB_DIR = os.getcwd()
OUT_DIR = os.path.join(os.getcwd(), 'SENT2_initiation_angles')

ODB_CASES = [
    ('GC', 4,  'SENT-II-GC-Miehe-R0-IPC4_UMAT_UMATHT-load.odb'),
    ('GC', 8,  'SENT-II-GC-Miehe-R0-IPC8_UMAT_UMATHT-load.odb'),
    ('GC', 16, 'SENT-II-GC-Miehe-R0-IPC16_UMAT_UMATHT-load.odb'),
    ('LE', 4,  'SENT-II-LE-Miehe-R0-IPC4_UMAT_UMATHT-load.odb'),
    ('LE', 8,  'SENT-II-LE-Miehe-R0-IPC8_UMAT_UMATHT-load.odb'),
    ('LE', 16, 'SENT-II-LE-Miehe-R0-IPC16_UMAT_UMATHT-load.odb'),
]

INSTANCE_NAME = 'PART-1-1'
STEP_NAME = 'Step-1'

# Geometric initial crack tip and initial crack direction.
# Current SENT-II model: horizontal initial crack ending at (0.5, 0.5).
INITIAL_TIP_X = 0.5
INITIAL_TIP_Y = 0.5
INITIAL_CRACK_ANGLE_DEG = 0.0

# 1 unit of Step time = 1 fatigue cycle
CYCLES_PER_TIME = 1.0

# To avoid reading all IPC frames, approximately one frame per cycle is sampled.
FRAME_SAMPLE_INTERVAL_CYCLES = 1.0

# Main phase-field threshold used to decide when a sufficiently long
# early crack has formed.
FRAME_SELECTION_THRESHOLD = 0.90

# The first sampled frame with crack centerline length >= this value is used.
# It only needs to be longer than ANGLE_FIT_LENGTH.
TARGET_EARLY_CRACK_LENGTH_MM = 0.060

# Local arc length used for fitting the initiation direction.
ANGLE_FIT_LENGTH_MM = 0.030

# Threshold sensitivity study. Main result can use 0.90;
# mean/std across these values can be used as a small error bar.
PHI_THRESHOLDS = [0.85, 0.90, 0.95]

# y-slice thickness for extracting a centerline from the finite-width
# phase-field crack band.
BIN_WIDTH_MM = 0.0015

# Light smoothing of the extracted centerline.
SMOOTH_WINDOW = 5

# Search window around the initial crack tip. This prevents remote high-phi
# regions from contaminating initiation-angle extraction.
SEARCH_X_MIN = 0.495
SEARCH_X_MAX = 0.620
SEARCH_Y_MIN = 0.360
SEARCH_Y_MAX = 0.505

# Pure Mode-II maximum tangential stress (MTS) prediction
MTS_MODE_II_ANGLE_DEG = 70.528779

# ======================================================================


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def scalar_value(data):
    try:
        return float(data)
    except (TypeError, ValueError):
        return float(data[0])


def get_phase_output(frame):
    if 'NT11' in frame.fieldOutputs:
        return frame.fieldOutputs['NT11'], 'NT11'
    if 'NT' in frame.fieldOutputs:
        return frame.fieldOutputs['NT'], 'NT'
    raise KeyError('Neither NT11 nor NT found at frameValue=%s' %
                   str(frame.frameValue))


def moving_average(values, window):
    n = len(values)
    if window < 3 or n < window:
        return list(values)

    half = window // 2
    padded = [values[0]] * half + list(values) + [values[-1]] * half
    out = []
    for i in range(n):
        out.append(sum(padded[i:i + window]) / float(window))
    return out


def extract_phi_nodes(frame, inst, coords):
    fld, field_name = get_phase_output(frame)
    vals = fld.getSubset(region=inst).values

    rows = []
    for v in vals:
        if not hasattr(v, 'nodeLabel'):
            continue

        xyz = coords.get(v.nodeLabel)
        if xyz is None:
            continue

        x, y, z = xyz

        if x < SEARCH_X_MIN or x > SEARCH_X_MAX:
            continue
        if y < SEARCH_Y_MIN or y > SEARCH_Y_MAX:
            continue

        phi = scalar_value(v.data)
        rows.append((v.nodeLabel, x, y, z, phi))

    return rows, field_name


def centerline_from_rows(rows, threshold):
    """
    Extract centerline using y-binning and phase-field-weighted centroids.
    The method assumes a single crack that runs monotonically downward, as in
    the shear-loaded SENT specimen.
    """
    selected = []
    for label, x, y, z, phi in rows:
        if phi >= threshold:
            selected.append((x, y, phi))

    if len(selected) < 4:
        return None

    y_min = min(p[1] for p in selected)

    bins = {}
    for x, y, phi in selected:
        bid = int(math.floor((y - y_min) / BIN_WIDTH_MM))
        if bid not in bins:
            bins[bid] = []
        bins[bid].append((x, y, phi))

    centers = []
    for bid in sorted(bins.keys(), reverse=True):
        pts = bins[bid]

        sw = 0.0
        sx = 0.0
        sy = 0.0

        for x, y, phi in pts:
            w = (phi - threshold) + 0.01
            sw += w
            sx += w * x
            sy += w * y

        if sw > 0.0:
            centers.append((sx / sw, sy / sw))

    if len(centers) < 4:
        return None

    # Sort from the initial tip towards the propagating crack.
    centers.sort(key=lambda p: p[1], reverse=True)

    # Anchor the geometric initial crack tip.
    xvals = [INITIAL_TIP_X]
    yvals = [INITIAL_TIP_Y]

    for x, y in centers:
        # Avoid an artificial duplicate exactly at the initial tip.
        if math.hypot(x - INITIAL_TIP_X, y - INITIAL_TIP_Y) > 1.0e-6:
            xvals.append(x)
            yvals.append(y)

    if len(xvals) < 4:
        return None

    xsm = moving_average(xvals, SMOOTH_WINDOW)
    ysm = moving_average(yvals, SMOOTH_WINDOW)

    # Keep the exact initial tip after smoothing.
    xsm[0] = INITIAL_TIP_X
    ysm[0] = INITIAL_TIP_Y

    # cumulative arc length
    s = [0.0]
    for i in range(1, len(xsm)):
        ds = math.hypot(xsm[i] - xsm[i - 1],
                        ysm[i] - ysm[i - 1])
        s.append(s[-1] + ds)

    return {
        'x_raw': xvals,
        'y_raw': yvals,
        'x': xsm,
        'y': ysm,
        's': s,
        'length': s[-1]
    }


def fit_initiation_angle(centerline):
    """
    Orthogonal/PCA line fit over 0 <= s <= ANGLE_FIT_LENGTH_MM.
    Returns signed angle relative to the +x direction.
    """
    x = []
    y = []

    for xi, yi, si in zip(centerline['x'],
                          centerline['y'],
                          centerline['s']):
        if si <= ANGLE_FIT_LENGTH_MM + 1.0e-12:
            x.append(xi)
            y.append(yi)

    if len(x) < 3:
        return None

    pts = np.array(list(zip(x, y)), dtype=float)
    ctr = pts.mean(axis=0)
    aa = pts - ctr

    # SVD / PCA direction
    u, sval, vh = np.linalg.svd(aa, full_matrices=False)
    direction = vh[0, :]

    # Orient along crack propagation
    chord = pts[-1, :] - pts[0, :]
    if np.dot(direction, chord) < 0.0:
        direction = -direction

    theta_global = math.degrees(
        math.atan2(direction[1], direction[0])
    )

    # Relative to original crack direction.
    theta_rel = theta_global - INITIAL_CRACK_ANGLE_DEG

    # Wrap into [-180, 180]
    while theta_rel > 180.0:
        theta_rel -= 360.0
    while theta_rel < -180.0:
        theta_rel += 360.0

    return theta_rel


def sampled_frames(step):
    """
    Return approximately one frame per requested cycle interval.
    """
    selected = []
    next_cycle = 0.0

    for i, frame in enumerate(step.frames):
        cycle = float(frame.frameValue) * CYCLES_PER_TIME

        if cycle + 1.0e-9 >= next_cycle:
            selected.append((i, frame, cycle))
            next_cycle += FRAME_SAMPLE_INTERVAL_CYCLES

    # Always retain the last frame
    if step.frames:
        last_i = len(step.frames) - 1
        last = step.frames[last_i]
        last_cycle = float(last.frameValue) * CYCLES_PER_TIME

        if len(selected) == 0 or selected[-1][0] != last_i:
            selected.append((last_i, last, last_cycle))

    return selected


def select_early_crack_frame(step, inst, coords):
    """
    Find the earliest sampled frame where the crack length extracted using
    FRAME_SELECTION_THRESHOLD reaches TARGET_EARLY_CRACK_LENGTH_MM.
    """
    frames = sampled_frames(step)

    print('  sampled frames to inspect:', len(frames))

    last_valid = None

    for count, (frame_index, frame, cycle) in enumerate(frames):
        rows, field_name = extract_phi_nodes(frame, inst, coords)
        cl = centerline_from_rows(rows, FRAME_SELECTION_THRESHOLD)

        if cl is None:
            continue

        last_valid = (frame_index, frame, cycle, rows, field_name, cl)

        if cl['length'] >= TARGET_EARLY_CRACK_LENGTH_MM:
            print('  selected frame:')
            print('    cycle        = %.6f' % cycle)
            print('    frame index  = %d' % frame_index)
            print('    crack length = %.6f mm' % cl['length'])
            return last_valid

    if last_valid is None:
        raise RuntimeError(
            'No valid phase-field crack centerline was found.'
        )

    print('  WARNING: target early crack length was not reached.')
    print('  Using last valid sampled frame:')
    print('    cycle        = %.6f' % last_valid[2])
    print('    frame index  = %d' % last_valid[0])
    print('    crack length = %.6f mm' % last_valid[5]['length'])

    return last_valid


def write_centerline_csv(path, model, ipc, cycle, frame_index,
                         threshold, centerline):
    mode = 'wb' if sys.version_info[0] < 3 else 'w'

    with open(path, mode) as f:
        w = csv.writer(f)
        w.writerow([
            'model', 'IPC', 'cycle_N', 'frame_index',
            'phi_threshold',
            'point_index',
            's_mm',
            'x_mm', 'y_mm',
            'dx_from_initial_tip_mm',
            'dy_from_initial_tip_mm'
        ])

        for i, (s, x, y) in enumerate(zip(centerline['s'],
                                          centerline['x'],
                                          centerline['y'])):
            w.writerow([
                model, ipc, cycle, frame_index,
                threshold,
                i, s, x, y,
                x - INITIAL_TIP_X,
                y - INITIAL_TIP_Y
            ])


def process_case(model, ipc, odb_name):
    odb_path = os.path.join(ODB_DIR, odb_name)

    if not os.path.isfile(odb_path):
        raise IOError('ODB not found: %s' % odb_path)

    print('')
    print('=' * 72)
    print('Processing:', odb_name)
    print('Model:', model, 'IPC:', ipc)
    print('=' * 72)

    odb = openOdb(path=odb_path, readOnly=True)

    try:
        if STEP_NAME not in odb.steps:
            raise KeyError('Step not found: %s' % STEP_NAME)

        if INSTANCE_NAME not in odb.rootAssembly.instances:
            names = sorted(odb.rootAssembly.instances.keys())
            raise KeyError(
                'Instance "%s" not found. Available: %s'
                % (INSTANCE_NAME, ', '.join(names))
            )

        step = odb.steps[STEP_NAME]
        inst = odb.rootAssembly.instances[INSTANCE_NAME]

        coords = {}
        for n in inst.nodes:
            xyz = n.coordinates
            x = float(xyz[0])
            y = float(xyz[1]) if len(xyz) > 1 else 0.0
            z = float(xyz[2]) if len(xyz) > 2 else 0.0
            coords[n.label] = (x, y, z)

        selected = select_early_crack_frame(step, inst, coords)

        frame_index, frame, cycle, rows, field_name, selection_cl = selected

        result = {
            'model': model,
            'IPC': ipc,
            'odb_name': odb_name,
            'selected_cycle_N': cycle,
            'selected_frame_index': frame_index,
            'selected_frame_value': float(frame.frameValue),
            'field_name': field_name,
            'selection_threshold': FRAME_SELECTION_THRESHOLD,
            'selection_crack_length_mm': selection_cl['length'],
        }

        angles = []

        for threshold in PHI_THRESHOLDS:
            cl = centerline_from_rows(rows, threshold)

            if cl is None:
                result['angle_%03d_deg' % int(round(threshold * 100))] = ''
                continue

            theta = fit_initiation_angle(cl)

            if theta is None:
                result['angle_%03d_deg' % int(round(threshold * 100))] = ''
                continue

            mag = abs(theta)
            angles.append(mag)

            key = 'angle_%03d_deg' % int(round(threshold * 100))
            result[key] = mag

            result['crack_length_%03d_mm' %
                   int(round(threshold * 100))] = cl['length']

            csv_name = '%s_IPC%d_phi%03d_centerline.csv' % (
                model, ipc, int(round(threshold * 100))
            )
            write_centerline_csv(
                os.path.join(OUT_DIR, csv_name),
                model, ipc, cycle, frame_index,
                threshold, cl
            )

            print('  phi >= %.2f: |theta0| = %.3f deg, L = %.5f mm' %
                  (threshold, mag, cl['length']))

        if len(angles) > 0:
            mean_angle = sum(angles) / float(len(angles))
            var = sum((a - mean_angle) ** 2 for a in angles) / float(len(angles))
            std_angle = math.sqrt(var)

            result['initiation_angle_mean_deg'] = mean_angle
            result['initiation_angle_std_deg'] = std_angle
            result['deviation_from_MTS_deg'] = abs(
                mean_angle - MTS_MODE_II_ANGLE_DEG
            )
        else:
            result['initiation_angle_mean_deg'] = ''
            result['initiation_angle_std_deg'] = ''
            result['deviation_from_MTS_deg'] = ''

        return result

    finally:
        odb.close()


def write_summary(rows):
    summary_path = os.path.join(
        OUT_DIR, 'SENT2_initiation_angle_summary.csv'
    )

    fields = [
        'model',
        'IPC',
        'odb_name',
        'selected_cycle_N',
        'selected_frame_index',
        'selected_frame_value',
        'field_name',
        'selection_threshold',
        'selection_crack_length_mm',
        'angle_085_deg',
        'angle_090_deg',
        'angle_095_deg',
        'initiation_angle_mean_deg',
        'initiation_angle_std_deg',
        'deviation_from_MTS_deg',
        'crack_length_085_mm',
        'crack_length_090_mm',
        'crack_length_095_mm',
    ]

    mode = 'wb' if sys.version_info[0] < 3 else 'w'

    with open(summary_path, mode) as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()

        for row in rows:
            clean = {}
            for field in fields:
                clean[field] = row.get(field, '')
            w.writerow(clean)

    # Origin-friendly compact table
    origin_path = os.path.join(
        OUT_DIR, 'SENT2_initiation_angle_Origin.csv'
    )

    with open(origin_path, mode) as f:
        w = csv.writer(f)
        w.writerow([
            'IPC',
            'GC_angle_mean_deg',
            'GC_angle_std_deg',
            'LE_angle_mean_deg',
            'LE_angle_std_deg',
            'MTS_theory_deg'
        ])

        ipcs = sorted(set(row['IPC'] for row in rows))

        for ipc in ipcs:
            gc = None
            le = None

            for row in rows:
                if row['IPC'] == ipc and row['model'] == 'GC':
                    gc = row
                if row['IPC'] == ipc and row['model'] == 'LE':
                    le = row

            w.writerow([
                ipc,
                gc.get('initiation_angle_mean_deg', '') if gc else '',
                gc.get('initiation_angle_std_deg', '') if gc else '',
                le.get('initiation_angle_mean_deg', '') if le else '',
                le.get('initiation_angle_std_deg', '') if le else '',
                MTS_MODE_II_ANGLE_DEG
            ])

    print('')
    print('Summary written:')
    print('  ', summary_path)
    print('Origin table:')
    print('  ', origin_path)


def main():
    ensure_dir(OUT_DIR)

    print('=' * 72)
    print('SENT-II initiation-angle batch extractor')
    print('ODB directory   :', ODB_DIR)
    print('Output directory:', OUT_DIR)
    print('Theory |theta0| : %.6f deg (pure Mode-II MTS)' %
          MTS_MODE_II_ANGLE_DEG)
    print('=' * 72)

    rows = []

    for model, ipc, odb_name in ODB_CASES:
        try:
            row = process_case(model, ipc, odb_name)
            rows.append(row)
        except Exception as exc:
            print('')
            print('ERROR processing:', odb_name)
            print(str(exc))

    if rows:
        write_summary(rows)

    print('')
    print('=' * 72)
    print('Done. Files in ' + OUT_DIR)
    print('=' * 72)


if __name__ == '__main__':
    main()
