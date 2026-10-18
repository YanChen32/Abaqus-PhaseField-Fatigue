# -*- coding: utf-8 -*-
"""
Nodal phase field of a GC and an LE result of the shear-loaded SENT specimen at matched crack extension
(article Fig. 6(a)).

The crack paths are compared at about the same crack extension, not at the same cycle number. The stages below
are read from the a(N) curves of the article at IPC = 16:
    a ~= 0.10 mm: GC N ~= 51,  LE N ~= 39
    a ~= 0.20 mm: GC N ~= 95,  LE N ~= 76
    a ~= 0.30 mm: GC N ~= 145, LE N ~= 119
    a ~= 0.40 mm: GC N ~= 213, LE N ~= 176
Other results need other stages (edit MATCHED_STAGES).

Usage:
    abaqus python export_sent2_path_raw.py GC_job.odb LE_job.odb

Output (folder SENT2_path_raw in the working directory):
    mesh_nodes.csv, mesh_elements.csv, manifest.csv, GC_a010_N0051.csv, LE_a010_N0039.csv, ...
"""

from __future__ import print_function
import os
import sys
import csv
from odbAccess import openOdb

# ----------------------------------------------------------------------
# User settings
# ----------------------------------------------------------------------

if len(sys.argv) != 3:
    raise SystemExit('Usage: abaqus python export_sent2_path_raw.py GC_job.odb LE_job.odb')
GC_ODB = os.path.abspath(sys.argv[1])
LE_ODB = os.path.abspath(sys.argv[2])

OUT_DIR = os.path.join(os.getcwd(), 'SENT2_path_raw')

INSTANCE_NAME = 'PART-1-1'
STEP_NAME = 'Step-1'

# Current model convention:
# 1 unit of Step time = 1 fatigue cycle
CYCLES_PER_TIME = 1.0

# Compare crack PATH at approximately equal crack extension.
# Each tuple: (label, target_a_mm, GC_cycle, LE_cycle)
MATCHED_STAGES = [
    ('a010', 0.10,  51.0,  39.0),
    ('a020', 0.20,  95.0,  76.0),
    ('a030', 0.30, 145.0, 119.0),
    ('a040', 0.40, 213.0, 176.0),
]

# Optional same-cycle snapshots. Useful for showing kinetics,
# but NOT for judging whether the crack trajectories are the same.
SAME_CYCLE_STAGES = []   # cycle numbers, e.g. [100.0, 150.0, 250.0]

# Export all nodal phase-field values in the selected frames.
# This avoids introducing a threshold during ODB extraction.
EXPORT_ALL_NODES = True

# Used only if EXPORT_ALL_NODES = False.
PHI_EXPORT_MIN = 0.02

# ----------------------------------------------------------------------


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def open_model(odb_path):
    if not os.path.isfile(odb_path):
        raise IOError('ODB not found: %s' % odb_path)

    odb = openOdb(path=odb_path, readOnly=True)

    if STEP_NAME not in odb.steps:
        odb.close()
        raise KeyError('Step "%s" not found in %s' % (STEP_NAME, odb_path))

    if INSTANCE_NAME not in odb.rootAssembly.instances:
        names = sorted(odb.rootAssembly.instances.keys())
        odb.close()
        raise KeyError(
            'Instance "%s" not found. Available instances: %s'
            % (INSTANCE_NAME, ', '.join(names))
        )

    step = odb.steps[STEP_NAME]
    inst = odb.rootAssembly.instances[INSTANCE_NAME]
    return odb, step, inst


def get_phase_output(frame):
    if 'NT11' in frame.fieldOutputs:
        return frame.fieldOutputs['NT11'], 'NT11'
    if 'NT' in frame.fieldOutputs:
        return frame.fieldOutputs['NT'], 'NT'
    raise KeyError('Neither NT11 nor NT exists in frame %s' % str(frame.frameValue))


def scalar_value(data):
    try:
        return float(data)
    except (TypeError, ValueError):
        return float(data[0])


def nearest_frame(step, target_cycle):
    target_time = float(target_cycle) / CYCLES_PER_TIME
    best_i = None
    best_frame = None
    best_err = None

    for i, frame in enumerate(step.frames):
        err = abs(float(frame.frameValue) - target_time)
        if best_err is None or err < best_err:
            best_err = err
            best_i = i
            best_frame = frame

    return best_i, best_frame, best_err


def export_mesh(inst, out_dir):
    nodes_path = os.path.join(out_dir, 'mesh_nodes.csv')
    elems_path = os.path.join(out_dir, 'mesh_elements.csv')

    mode = 'wb' if sys.version_info[0] < 3 else 'w'

    with open(nodes_path, mode) as f:
        w = csv.writer(f)
        w.writerow(['node_label', 'x_mm', 'y_mm', 'z_mm'])
        for n in inst.nodes:
            xyz = n.coordinates
            x = float(xyz[0])
            y = float(xyz[1]) if len(xyz) > 1 else 0.0
            z = float(xyz[2]) if len(xyz) > 2 else 0.0
            w.writerow([n.label, x, y, z])

    max_conn = 0
    for e in inst.elements:
        if len(e.connectivity) > max_conn:
            max_conn = len(e.connectivity)

    with open(elems_path, mode) as f:
        w = csv.writer(f)
        header = ['element_label'] + [
            'node_%d' % (i + 1) for i in range(max_conn)
        ]
        w.writerow(header)

        for e in inst.elements:
            row = [e.label] + list(e.connectivity)
            if len(e.connectivity) < max_conn:
                row += [''] * (max_conn - len(e.connectivity))
            w.writerow(row)

    print('Wrote mesh:')
    print('  ' + nodes_path)
    print('  ' + elems_path)


def export_frame(model_name, stage_label, target_a, target_cycle,
                 step, inst, out_dir):
    frame_index, frame, frame_err = nearest_frame(step, target_cycle)

    phase, phase_name = get_phase_output(frame)
    vals = phase.getSubset(region=inst).values

    coords = {}
    for n in inst.nodes:
        xyz = n.coordinates
        coords[n.label] = (
            float(xyz[0]),
            float(xyz[1]) if len(xyz) > 1 else 0.0,
            float(xyz[2]) if len(xyz) > 2 else 0.0
        )

    actual_cycle = float(frame.frameValue) * CYCLES_PER_TIME

    filename = '%s_%s_N%04d.csv' % (
        model_name,
        stage_label,
        int(round(actual_cycle))
    )
    out_path = os.path.join(out_dir, filename)

    mode = 'wb' if sys.version_info[0] < 3 else 'w'

    n_written = 0
    with open(out_path, mode) as f:
        w = csv.writer(f)
        w.writerow([
            'model',
            'stage_label',
            'target_crack_extension_mm',
            'target_cycle',
            'actual_cycle',
            'frame_index',
            'frame_value',
            'field_name',
            'node_label',
            'x_mm',
            'y_mm',
            'z_mm',
            'phi'
        ])

        for v in vals:
            if not hasattr(v, 'nodeLabel'):
                continue

            phi = scalar_value(v.data)

            if (not EXPORT_ALL_NODES) and phi < PHI_EXPORT_MIN:
                continue

            xyz = coords.get(v.nodeLabel)
            if xyz is None:
                continue

            w.writerow([
                model_name,
                stage_label,
                target_a,
                target_cycle,
                actual_cycle,
                frame_index,
                float(frame.frameValue),
                phase_name,
                v.nodeLabel,
                xyz[0],
                xyz[1],
                xyz[2],
                phi
            ])
            n_written += 1

    print('%s %s' % (model_name, stage_label))
    print('  target a       : %.6f mm' % float(target_a))
    print('  target cycle   : %.6f' % float(target_cycle))
    print('  actual cycle   : %.6f' % actual_cycle)
    print('  frame index    : %d' % frame_index)
    print('  frame mismatch : %.6e step-time' % frame_err)
    print('  nodes written  : %d' % n_written)
    print('  file           : %s' % out_path)

    return {
        'model': model_name,
        'stage_label': stage_label,
        'target_a_mm': target_a,
        'target_cycle': target_cycle,
        'actual_cycle': actual_cycle,
        'frame_index': frame_index,
        'frame_value': float(frame.frameValue),
        'frame_time_error': frame_err,
        'field_name': phase_name,
        'nodes_written': n_written,
        'csv_file': filename
    }


def export_same_cycle(model_name, target_cycle, step, inst, out_dir):
    stage_label = 'sameN%04d' % int(round(target_cycle))
    return export_frame(
        model_name,
        stage_label,
        0.0,
        target_cycle,
        step,
        inst,
        out_dir
    )


def write_manifest(rows, out_dir):
    path = os.path.join(out_dir, 'manifest.csv')
    mode = 'wb' if sys.version_info[0] < 3 else 'w'

    fields = [
        'model',
        'stage_label',
        'target_a_mm',
        'target_cycle',
        'actual_cycle',
        'frame_index',
        'frame_value',
        'frame_time_error',
        'field_name',
        'nodes_written',
        'csv_file'
    ]

    with open(path, mode) as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    print('Wrote manifest:')
    print('  ' + path)


def main():
    ensure_dir(OUT_DIR)

    print('=' * 72)
    print('SENT-II GC/LE raw crack-path exporter')
    print('Output directory: ' + OUT_DIR)
    print('=' * 72)

    gc_odb, gc_step, gc_inst = open_model(GC_ODB)
    le_odb, le_step, le_inst = open_model(LE_ODB)

    try:
        print('GC nodes/elements: %d %d' % (len(gc_inst.nodes), len(gc_inst.elements)))
        print('LE nodes/elements: %d %d' % (len(le_inst.nodes), len(le_inst.elements)))

        if len(gc_inst.nodes) != len(le_inst.nodes):
            print('WARNING: GC and LE node counts differ.')
        if len(gc_inst.elements) != len(le_inst.elements):
            print('WARNING: GC and LE element counts differ.')

        # Controlled comparison should use the same mesh.
        # Export the GC mesh once for later graph/centerline analysis.
        export_mesh(gc_inst, OUT_DIR)

        manifest = []

        for stage_label, target_a, gc_cycle, le_cycle in MATCHED_STAGES:
            manifest.append(
                export_frame(
                    'GC', stage_label, target_a, gc_cycle,
                    gc_step, gc_inst, OUT_DIR
                )
            )
            manifest.append(
                export_frame(
                    'LE', stage_label, target_a, le_cycle,
                    le_step, le_inst, OUT_DIR
                )
            )

        for target_cycle in SAME_CYCLE_STAGES:
            manifest.append(
                export_same_cycle(
                    'GC', target_cycle,
                    gc_step, gc_inst, OUT_DIR
                )
            )
            manifest.append(
                export_same_cycle(
                    'LE', target_cycle,
                    le_step, le_inst, OUT_DIR
                )
            )

        write_manifest(manifest, OUT_DIR)

    finally:
        gc_odb.close()
        le_odb.close()

    print('=' * 72)
    print('Done. Files in ' + OUT_DIR)
    print('=' * 72)


if __name__ == '__main__':
    main()
