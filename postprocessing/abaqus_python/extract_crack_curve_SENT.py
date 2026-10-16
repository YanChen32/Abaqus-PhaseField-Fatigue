# -*- coding: utf-8 -*-
"""
Crack extension a(N) of the Mode-I SENT specimen (article Fig. 2).

The crack front is the largest x coordinate of the nodes with phi >= 0.95; a is measured from the initial
crack tip at x = 0.5 mm. The CSV is written next to the ODB name in the working directory.

Usage:
    abaqus python extract_crack_curve_SENT.py jobname.odb
"""
from __future__ import print_function
import os, sys, csv
from odbAccess import openOdb

PHI_THRESHOLD = 0.95
INITIAL_CRACK_TIP_X = 0.5
INSTANCE_NAME = 'PART-1-1'
STEP_NAME = 'Step-1'

if len(sys.argv) < 2:
    raise SystemExit('Usage: abaqus python extract_crack_curve_SENT.py jobname.odb')

odb_path = sys.argv[1]
odb = openOdb(path=odb_path, readOnly=True)
step = odb.steps[STEP_NAME]
inst = odb.rootAssembly.instances[INSTANCE_NAME]
coords = dict((n.label, n.coordinates) for n in inst.nodes)

rows = []
for frame in step.frames:
    if 'NT11' in frame.fieldOutputs:
        fld = frame.fieldOutputs['NT11']
    elif 'NT' in frame.fieldOutputs:
        fld = frame.fieldOutputs['NT']
    else:
        continue

    vals = fld.getSubset(region=inst).values
    xmax = None
    for v in vals:
        try:
            phi = float(v.data)
        except TypeError:
            phi = float(v.data[0])

        if phi >= PHI_THRESHOLD:
            xyz = coords.get(v.nodeLabel)
            if xyz is not None:
                x = float(xyz[0])
                if x >= INITIAL_CRACK_TIP_X - 1.0e-8:
                    if xmax is None or x > xmax:
                        xmax = x

    if xmax is None:
        a = 0.0
    else:
        a = max(0.0, xmax - INITIAL_CRACK_TIP_X)

    N = float(frame.frameValue)
    rows.append((N, a))

odb.close()
base = os.path.splitext(os.path.basename(odb_path))[0]
out_csv = base + '_a_vs_N.csv'
mode = 'wb' if sys.version_info[0] < 3 else 'w'
with open(out_csv, mode) as f:
    w = csv.writer(f)
    w.writerow(['N_cycles', 'crack_extension_a_mm'])
    for row in rows:
        w.writerow(row)

print('Wrote:', out_csv)
