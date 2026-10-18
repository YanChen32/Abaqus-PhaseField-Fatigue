# -*- coding: utf-8 -*-
"""
Crack-tip history of a half compact-tension (CT) model computed with the cycle-domain subroutine.

For every ODB in the working directory (optionally filtered by ODB_NAME_FILTER) the script
1. reads NT11 on the crack line (symmetry line y = CRACK_LINE_Y);
2. takes the crack tip as the crossing NT11 = 0.95, kept monotonic in time;
3. writes <job>_crack_history.csv with the columns frame, N, x_tip_raw_mm, x_tip_monotonic_mm, delta_a_mm, ...;
4. writes CT_crack_history_summary.csv with the mean growth rate between 5 l0 and 15 l0 of crack extension.

In the cycle-domain subroutine the step time is the cycle number, so frame.frameValue is used as N.
The crack-history files are the input of postprocessing/identification/paris_identify.py.

Usage:
    abaqus python extract_CT_crack_history.py
or in Abaqus/CAE: File -> Run Script.
"""

from __future__ import print_function

import os
import re
import math
import traceback

from odbAccess import openOdb


# ======================================================================
# USER SETTINGS -- normally edit only this block
# ======================================================================

# Folder containing the ODB files (default: working directory)
ODB_DIR = os.getcwd()

# Process only ODB files whose filename contains this text.
# '' processes every .odb in ODB_DIR.
ODB_NAME_FILTER = ''

# Step to read.
# Use '' to automatically use the last analysis step.
STEP_NAME = 'Step-1'

# Instance name.
# Use '' to automatically choose the instance with the largest number of nodes.
INSTANCE_NAME = ''

# Phase-field output variable
FIELD_NAME = 'NT11'

# Crack definition:
# crack is considered fully developed behind the point where NT11 crosses this.
PHI_THRESHOLD = 0.95

# y coordinate of the crack line (symmetry line of the half model).
CRACK_LINE_Y = 0.0

# Nodes within this distance from CRACK_LINE_Y are used.
# With l0 = 0.1 mm, 0.02 mm is a sensible starting value.
# If the script reports too few crack-line nodes, increase this slightly.
Y_TOL = 0.02

# Optional X search window for the crack path.
# Use None to impose no limit.
X_MIN_SEARCH = None
X_MAX_SEARCH = None

# Direction of crack growth: +1 = increasing x, -1 = decreasing x.
CRACK_GROWTH_DIRECTION = +1

# Phase-field length scale l0, mm
L0 = 0.1

# Crack-extension window (in units of l0) of the mean rate in the summary file
A1_FACTOR = 5.0
A2_FACTOR = 15.0

# Read every frame. Change to 2, 5, ... only if the ODB has extremely many frames.
FRAME_STRIDE = 1

# Additional cycle offset if the step starts after another fatigue step.
CYCLE_OFFSET = 0.0

# Whether to enforce irreversible crack-tip motion in post-processing.
# Recommended because tiny NT11 numerical fluctuations can otherwise make
# the detected crack tip move backwards by one node.
ENFORCE_MONOTONIC_TIP = True


# ======================================================================
# HELPER FUNCTIONS
# ======================================================================

def scalar_value(data):
    """Convert an Abaqus scalar FieldValue.data to float."""
    try:
        return float(data)
    except Exception:
        try:
            return float(data[0])
        except Exception:
            return float('nan')


def find_load_from_name(filename):
    """Extract LOAD_0400 -> 400 from filename, if available."""
    m = re.search(r'LOAD[_\-]?(\d+)', filename.upper())
    if m:
        try:
            return int(m.group(1))
        except Exception:
            pass
    return None


def choose_instance(odb):
    """Choose requested instance or automatically choose the largest instance."""
    root = odb.rootAssembly

    if INSTANCE_NAME:
        if INSTANCE_NAME not in root.instances.keys():
            raise RuntimeError(
                'Instance "%s" not found. Available instances: %s'
                % (INSTANCE_NAME, str(root.instances.keys()))
            )
        return root.instances[INSTANCE_NAME]

    best = None
    best_n = -1
    for name in root.instances.keys():
        inst = root.instances[name]
        n = len(inst.nodes)
        if n > best_n:
            best_n = n
            best = inst

    if best is None:
        raise RuntimeError('No instance with nodes was found in the ODB.')

    return best


def build_node_coordinate_map(inst):
    """Map nodeLabel -> (x, y, z)."""
    out = {}
    for node in inst.nodes:
        c = node.coordinates
        if len(c) == 2:
            out[node.label] = (float(c[0]), float(c[1]), 0.0)
        else:
            out[node.label] = (float(c[0]), float(c[1]), float(c[2]))
    return out


def on_crack_line(x, y):
    if abs(y - CRACK_LINE_Y) > Y_TOL:
        return False
    if X_MIN_SEARCH is not None and x < X_MIN_SEARCH:
        return False
    if X_MAX_SEARCH is not None and x > X_MAX_SEARCH:
        return False
    return True


def interpolate_tip(points, threshold, direction):
    """
    Abaqus/CAE-safe version:
    Do NOT use generator expressions inside max/min/sum, because Abaqus/CAE
    may shadow Python built-ins with Abaqus functions that reject generators.

    points: list of (x, phi), already restricted to crack line.

    Returns:
        x_tip, n_above, phi_max
    """
    if not points:
        return float('nan'), 0, float('nan')

    points = sorted(points, key=lambda z: z[0])

    if direction < 0:
        points = list(reversed(points))

    n_above = 0
    phi_max = None

    for item in points:
        phi = item[1]

        if phi >= threshold:
            n_above += 1

        if phi_max is None or phi > phi_max:
            phi_max = phi

    crossings = []

    for i in range(len(points) - 1):
        x1 = points[i][0]
        p1 = points[i][1]
        x2 = points[i + 1][0]
        p2 = points[i + 1][1]

        # Moving in growth direction: cracked -> intact
        if p1 >= threshold and p2 < threshold:
            if abs(p2 - p1) < 1.0e-15:
                xc = x1
            else:
                t = (threshold - p1) / (p2 - p1)
                xc = x1 + t * (x2 - x1)
            crossings.append(xc)

    if crossings:
        x_tip = crossings[0]
        for xc in crossings[1:]:
            if direction > 0:
                if xc > x_tip:
                    x_tip = xc
            else:
                if xc < x_tip:
                    x_tip = xc
        return x_tip, n_above, phi_max

    # Fallback if there is no clean crossing
    cracked_x = []
    for item in points:
        if item[1] >= threshold:
            cracked_x.append(item[0])

    if cracked_x:
        x_tip = cracked_x[0]
        for xc in cracked_x[1:]:
            if direction > 0:
                if xc > x_tip:
                    x_tip = xc
            else:
                if xc < x_tip:
                    x_tip = xc
        return x_tip, n_above, phi_max

    return float('nan'), 0, phi_max

def interpolate_cycle_at_extension(rows, target_da):
    """
    rows contain:
      N, x_tip_raw, x_tip_mono, delta_a, ...
    Uses piecewise-linear interpolation to obtain N(delta_a = target_da).
    """
    valid = []
    for row in rows:
        N = row[0]
        da = row[3]
        if not (math.isnan(N) or math.isnan(da)):
            valid.append((N, da))

    if not valid:
        return float('nan')

    if target_da <= valid[0][1]:
        return valid[0][0]

    for i in range(len(valid) - 1):
        n1, a1 = valid[i]
        n2, a2 = valid[i + 1]

        if a1 <= target_da <= a2:
            if abs(a2 - a1) < 1.0e-15:
                return n2
            t = (target_da - a1) / (a2 - a1)
            return n1 + t * (n2 - n1)

    return float('nan')


def fmt(x):
    if x is None:
        return ''
    try:
        if math.isnan(x):
            return ''
    except Exception:
        pass
    return '%.12g' % x


# ======================================================================
# PROCESS ONE ODB
# ======================================================================

def process_odb(odb_path):
    print('')
    print('------------------------------------------------------------')
    print('Processing: %s' % odb_path)

    odb = openOdb(path=odb_path, readOnly=True)

    try:
        # Step selection
        if STEP_NAME:
            if STEP_NAME not in odb.steps.keys():
                raise RuntimeError(
                    'Step "%s" not found. Available steps: %s'
                    % (STEP_NAME, str(odb.steps.keys()))
                )
            step = odb.steps[STEP_NAME]
            used_step_name = STEP_NAME
        else:
            step_names = list(odb.steps.keys())
            if not step_names:
                raise RuntimeError('No analysis step found in ODB.')
            used_step_name = step_names[-1]
            step = odb.steps[used_step_name]

        inst = choose_instance(odb)
        coords = build_node_coordinate_map(inst)

        # Crack-line node count from geometry
        crack_line_labels = []
        for label in coords:
            x, y, z = coords[label]
            if on_crack_line(x, y):
                crack_line_labels.append(label)

        print('Step            : %s' % used_step_name)
        print('Instance        : %s' % inst.name)
        print('Crack-line nodes: %d' % len(crack_line_labels))
        print('NT threshold    : %.4f' % PHI_THRESHOLD)
        print('Crack line Y    : %.6g +/- %.6g' % (CRACK_LINE_Y, Y_TOL))

        if len(crack_line_labels) < 3:
            raise RuntimeError(
                'Too few nodes were found on the crack line. '
                'Check CRACK_LINE_Y and Y_TOL.'
            )

        rows = []
        previous_tip = None

        frame_indices = list(range(0, len(step.frames), FRAME_STRIDE))
        if frame_indices and frame_indices[-1] != len(step.frames) - 1:
            frame_indices.append(len(step.frames) - 1)

        for iframe in frame_indices:
            frame = step.frames[iframe]

            if FIELD_NAME not in frame.fieldOutputs.keys():
                print('WARNING: %s missing in frame %d; skipped.'
                      % (FIELD_NAME, iframe))
                continue

            field = frame.fieldOutputs[FIELD_NAME]

            # Restrict to selected instance if possible.
            try:
                field_inst = field.getSubset(region=inst)
                values = field_inst.values
            except Exception:
                values = field.values

            points = []

            for value in values:
                try:
                    label = value.nodeLabel
                except Exception:
                    continue

                if label not in coords:
                    continue

                x, y, z = coords[label]

                if not on_crack_line(x, y):
                    continue

                phi = scalar_value(value.data)
                if math.isnan(phi):
                    continue

                points.append((x, phi))

            x_tip_raw, n_above, phi_max = interpolate_tip(
                points, PHI_THRESHOLD, CRACK_GROWTH_DIRECTION
            )

            x_tip = x_tip_raw

            if ENFORCE_MONOTONIC_TIP and previous_tip is not None:
                if not math.isnan(x_tip):
                    if CRACK_GROWTH_DIRECTION > 0:
                        if x_tip < previous_tip:
                            x_tip = previous_tip
                    else:
                        if x_tip > previous_tip:
                            x_tip = previous_tip

            if not math.isnan(x_tip):
                previous_tip = x_tip

            N = float(frame.frameValue) + CYCLE_OFFSET

            # delta_a is added after initial tip is known
            rows.append([
                N,              # 0
                x_tip_raw,      # 1
                x_tip,          # 2
                float('nan'),   # 3 delta_a
                n_above,        # 4
                phi_max,        # 5
                iframe          # 6
            ])

        if not rows:
            raise RuntimeError('No valid frame data were extracted.')

        # Initial crack tip = first valid monotonic tip
        initial_tip = float('nan')
        for row in rows:
            if not math.isnan(row[2]):
                initial_tip = row[2]
                break

        if math.isnan(initial_tip):
            raise RuntimeError(
                'Could not detect an initial crack tip. '
                'Check FIELD_NAME, threshold, crack line Y, and Y_TOL.'
            )

        for row in rows:
            xt = row[2]
            if math.isnan(xt):
                row[3] = float('nan')
            else:
                if CRACK_GROWTH_DIRECTION > 0:
                    row[3] = xt - initial_tip
                else:
                    row[3] = initial_tip - xt

        # Stable-window targets
        da1 = A1_FACTOR * L0
        da2 = A2_FACTOR * L0

        N1 = interpolate_cycle_at_extension(rows, da1)
        N2 = interpolate_cycle_at_extension(rows, da2)

        if (not math.isnan(N1)) and (not math.isnan(N2)) and N2 > N1:
            v_stable = (da2 - da1) / (N2 - N1)
        else:
            v_stable = float('nan')

        # Final valid state
        final_row = None
        for row in reversed(rows):
            if not math.isnan(row[3]):
                final_row = row
                break

        final_N = final_row[0]
        final_tip = final_row[2]
        final_da = final_row[3]

        first_N = rows[0][0]

        if final_N > first_N:
            v_overall = final_da / (final_N - first_N)
        else:
            v_overall = float('nan')

        # Per-ODB CSV
        base = os.path.splitext(os.path.basename(odb_path))[0]
        out_csv = os.path.join(ODB_DIR, base + '_crack_history.csv')

        f = open(out_csv, 'w')
        try:
            f.write(
                'frame,N,x_tip_raw_mm,x_tip_monotonic_mm,'
                'delta_a_mm,n_nodes_phi_ge_threshold,max_NT11\n'
            )
            for row in rows:
                f.write(
                    '%d,%s,%s,%s,%s,%d,%s\n'
                    % (
                        row[6],
                        fmt(row[0]),
                        fmt(row[1]),
                        fmt(row[2]),
                        fmt(row[3]),
                        row[4],
                        fmt(row[5])
                    )
                )
        finally:
            f.close()

        print('Initial tip     : %.6g mm' % initial_tip)
        print('Final tip       : %.6g mm' % final_tip)
        print('Final delta_a   : %.6g mm' % final_da)
        print('N(5*l0)         : %s' % fmt(N1))
        print('N(15*l0)        : %s' % fmt(N2))
        print('v[5l0->15l0]    : %s mm/cycle' % fmt(v_stable))
        print('CSV written     : %s' % out_csv)

        return {
            'filename': os.path.basename(odb_path),
            'load_N': find_load_from_name(os.path.basename(odb_path)),
            'step': used_step_name,
            'instance': inst.name,
            'initial_tip': initial_tip,
            'final_N': final_N,
            'final_tip': final_tip,
            'final_da': final_da,
            'N_5l0': N1,
            'N_15l0': N2,
            'v_5l0_15l0': v_stable,
            'v_overall': v_overall,
            'history_csv': out_csv
        }

    finally:
        odb.close()


# ======================================================================
# MAIN
# ======================================================================

def main():
    if not os.path.isdir(ODB_DIR):
        raise RuntimeError('ODB_DIR does not exist: %s' % ODB_DIR)

    odb_files = []

    for name in os.listdir(ODB_DIR):
        low = name.lower()

        if not low.endswith('.odb'):
            continue

        if ODB_NAME_FILTER and ODB_NAME_FILTER.upper() not in name.upper():
            continue

        odb_files.append(name)

    if not odb_files:
        raise RuntimeError(
            'No matching ODB files found in:\n%s\nFilter: %s'
            % (ODB_DIR, ODB_NAME_FILTER)
        )

    # Sort by extracted load first, then filename
    def sort_key(name):
        load = find_load_from_name(name)
        if load is None:
            load = 10**12
        return (load, name)

    odb_files.sort(key=sort_key)

    print('============================================================')
    print('Half-CT phase-field crack-history extraction')
    print('============================================================')
    print('ODB directory  : %s' % ODB_DIR)
    print('Files found    : %d' % len(odb_files))
    print('Field          : %s' % FIELD_NAME)
    print('Threshold      : %.4f' % PHI_THRESHOLD)
    print('l0             : %.6g mm' % L0)
    print('Stable window  : %.3g*l0 -> %.3g*l0'
          % (A1_FACTOR, A2_FACTOR))
    print('============================================================')

    summaries = []

    for name in odb_files:
        path = os.path.join(ODB_DIR, name)
        try:
            result = process_odb(path)
            summaries.append(result)
        except Exception as err:
            print('')
            print('ERROR while processing %s' % name)
            print(str(err))
            print('--- Python traceback ---')
            traceback.print_exc()
            print('------------------------')
            print('This ODB was skipped.')

    # Summary CSV
    summary_path = os.path.join(ODB_DIR, 'CT_crack_history_summary.csv')

    f = open(summary_path, 'w')
    try:
        f.write(
            'odb,load_N,step,instance,initial_tip_mm,final_N,'
            'final_tip_mm,delta_a_final_mm,N_at_5l0,N_at_15l0,'
            'v_5l0_to_15l0_mm_per_cycle,v_overall_mm_per_cycle\n'
        )

        for s in summaries:
            f.write(
                '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'
                % (
                    s['filename'],
                    '' if s['load_N'] is None else str(s['load_N']),
                    s['step'],
                    s['instance'],
                    fmt(s['initial_tip']),
                    fmt(s['final_N']),
                    fmt(s['final_tip']),
                    fmt(s['final_da']),
                    fmt(s['N_5l0']),
                    fmt(s['N_15l0']),
                    fmt(s['v_5l0_15l0']),
                    fmt(s['v_overall'])
                )
            )
    finally:
        f.close()

    print('')
    print('============================================================')
    print('Finished.')
    print('Summary: %s' % summary_path)
    print('============================================================')


if __name__ == '__main__':
    main()
