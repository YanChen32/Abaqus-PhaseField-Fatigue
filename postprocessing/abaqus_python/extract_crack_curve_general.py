# -*- coding: utf-8 -*-
"""
Crack extension along a curved crack (shear-loaded SENT specimen, article Fig. 4).

The crack extension is the shortest-path distance, through the nodes with phi >= 0.90, from the initial crack
tip to the farthest node with phi >= 0.95. It is measured along the crack, so a kinked or curved crack is not
projected onto the x axis.

Usage
-----
1) One ODB by file name (looked up in ODB_DIR) or by full path:
   abaqus python extract_crack_curve_general.py Job-1.odb

2) All ODB files in ODB_DIR:
   abaqus python extract_crack_curve_general.py --all

ODB_DIR and RESULT_DIR default to the working directory.

Notes
-----
- The phase field is assumed to be stored as NT11 (or NT).
- frame.frameValue is assumed to be pseudo-time; CYCLES_PER_TIME converts it
  to cycle number. For the cycle-by-cycle SENT model, use 1.0.
- INITIAL_CRACK_TIP must be set to the physical initial crack-tip coordinate.
- The current implementation is intended for 2D continuum meshes. Connecting
  all nodes belonging to the same element gives diagonal graph links and
  avoids the strong Manhattan-distance error of edge-only graph tracking.
"""
from __future__ import print_function

import os
import sys
import csv
import math
import heapq
from odbAccess import openOdb


# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
ODB_DIR = os.getcwd()
RESULT_DIR = os.getcwd()

INSTANCE_NAME = 'PART-1-1'
STEP_NAME = 'Step-1'

# Crack definition. Use a high value for the crack core and a slightly lower
# value to keep the curved crack band connected.
TIP_THRESHOLD = 0.95
CONNECT_THRESHOLD = 0.90

# Modify this if the initial crack tip is not (0.5, 0.5).
INITIAL_CRACK_TIP = (0.5, 0.5)

# One unit of step time = one cycle for the present SENT cycle-by-cycle model.
CYCLES_PER_TIME = 1.0

# Extract once per this many cycles. Set to 0 or None to process every ODB frame.
SAMPLE_INTERVAL_N = 1.0

# Initial crack nodes are removed from the crack-growth candidates. A small
# neighborhood around the initial tip is kept so that the new crack remains
# connected to the original tip. Radius = SEED_RADIUS_FACTOR * local mesh size.
SEED_RADIUS_FACTOR = 3.0

# Numerical tolerance for coordinate and frame comparisons.
TOL = 1.0e-10


# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------
def _phi_scalar(data):
    try:
        return float(data)
    except (TypeError, ValueError):
        return float(data[0])


def _distance(a, b):
    n = min(len(a), len(b))
    s = 0.0
    for i in range(n):
        d = float(a[i]) - float(b[i])
        s += d * d
    return math.sqrt(s)


def _get_phase_field(frame, inst):
    if 'NT11' in frame.fieldOutputs:
        fld = frame.fieldOutputs['NT11']
    elif 'NT' in frame.fieldOutputs:
        fld = frame.fieldOutputs['NT']
    else:
        return {}

    result = {}
    vals = fld.getSubset(region=inst).values
    for v in vals:
        if not hasattr(v, 'nodeLabel'):
            continue
        result[v.nodeLabel] = _phi_scalar(v.data)
    return result


def _build_graph(inst, coords):
    """Build a weighted node graph from element connectivity.

    All node pairs in the same element are connected. For a 4-node quad this
    gives the usual 8-neighbor type connectivity, so a diagonal crack does not
    acquire the large Manhattan-length bias that occurs with edge-only links.
    """
    graph = dict((label, {}) for label in coords)

    for elem in inst.elements:
        conn = list(elem.connectivity)
        m = len(conn)
        for i in range(m):
            ni = conn[i]
            if ni not in coords:
                continue
            for j in range(i + 1, m):
                nj = conn[j]
                if nj not in coords:
                    continue
                w = _distance(coords[ni], coords[nj])
                if w <= TOL:
                    continue

                old = graph[ni].get(nj)
                if old is None or w < old:
                    graph[ni][nj] = w
                    graph[nj][ni] = w

    return graph


def _nearest_node(labels, coords, point):
    best = None
    best_d = None
    for label in labels:
        xyz = coords.get(label)
        if xyz is None:
            continue
        d = _distance(xyz, point)
        if best_d is None or d < best_d:
            best = label
            best_d = d
    return best, best_d


def _estimate_local_mesh_size(root, graph):
    if root is None or root not in graph:
        return 0.0
    lengths = list(graph[root].values())
    if not lengths:
        return 0.0
    lengths.sort()
    # Use the smallest local connection length; diagonal element links are
    # longer than the true edge size and should not define h.
    return float(lengths[0])


def _dijkstra(root, allowed, graph):
    if root is None or root not in allowed:
        return {}

    dist = {root: 0.0}
    heap = [(0.0, root)]

    while heap:
        du, u = heapq.heappop(heap)
        if du > dist.get(u, 1.0e300) + TOL:
            continue

        for v, w in graph.get(u, {}).items():
            if v not in allowed:
                continue
            nd = du + w
            if nd + TOL < dist.get(v, 1.0e300):
                dist[v] = nd
                heapq.heappush(heap, (nd, v))

    return dist


def _select_frames(frames, sample_interval_time):
    if sample_interval_time is None or sample_interval_time <= 0.0:
        return list(frames)

    if not frames:
        return []

    start = float(frames[0].frameValue)
    end = float(frames[-1].frameValue)
    selected = []
    used = set()

    target = start
    j = 0
    n = len(frames)

    while target <= end + TOL:
        while (j + 1 < n and
               abs(float(frames[j + 1].frameValue) - target) <=
               abs(float(frames[j].frameValue) - target)):
            j += 1

        if j not in used:
            selected.append(frames[j])
            used.add(j)
        target += sample_interval_time

    if (n - 1) not in used:
        selected.append(frames[-1])

    return selected


def _resolve_odb_path(arg):
    if os.path.isabs(arg) and os.path.isfile(arg):
        return arg

    if os.path.isfile(arg):
        return os.path.abspath(arg)

    name = arg
    if not name.lower().endswith('.odb'):
        name += '.odb'
    path = os.path.join(ODB_DIR, name)
    if os.path.isfile(path):
        return path

    raise IOError('ODB file not found: %s' % arg)


def _ensure_result_dir():
    if not os.path.isdir(RESULT_DIR):
        os.makedirs(RESULT_DIR)


# -----------------------------------------------------------------------------
# Crack extraction
# -----------------------------------------------------------------------------
def extract_one_odb(odb_path):
    print('Opening:', odb_path)
    odb = openOdb(path=odb_path, readOnly=True)

    try:
        if STEP_NAME not in odb.steps:
            raise KeyError('Step not found: %s' % STEP_NAME)
        step = odb.steps[STEP_NAME]

        if INSTANCE_NAME not in odb.rootAssembly.instances:
            raise KeyError('Instance not found: %s' % INSTANCE_NAME)
        inst = odb.rootAssembly.instances[INSTANCE_NAME]

        coords = dict((n.label, tuple(n.coordinates)) for n in inst.nodes)
        graph = _build_graph(inst, coords)

        if not step.frames:
            raise RuntimeError('No frames in step: %s' % STEP_NAME)

        # Initial crack state from the first frame.
        phi0 = _get_phase_field(step.frames[0], inst)
        if not phi0:
            raise RuntimeError('NT11/NT field not found in the first frame.')

        initial_core = set(label for label, phi in phi0.items()
                           if phi >= TIP_THRESHOLD)
        initial_connect = set(label for label, phi in phi0.items()
                              if phi >= CONNECT_THRESHOLD)

        # Root is the node nearest to the prescribed physical initial crack tip.
        search_labels = initial_connect if initial_connect else set(coords.keys())
        root, root_dist = _nearest_node(search_labels, coords, INITIAL_CRACK_TIP)
        if root is None:
            raise RuntimeError('Cannot find a node near INITIAL_CRACK_TIP.')

        h_local = _estimate_local_mesh_size(root, graph)
        if h_local <= 0.0:
            seed_radius = 0.0
        else:
            seed_radius = SEED_RADIUS_FACTOR * h_local

        seed_band = set()
        for label in initial_connect:
            if _distance(coords[label], INITIAL_CRACK_TIP) <= seed_radius + TOL:
                seed_band.add(label)
        seed_band.add(root)

        print('Initial crack tip target :', INITIAL_CRACK_TIP)
        print('Nearest mesh node        :', root, coords[root])
        print('Distance to target       :', root_dist)
        print('Local mesh size h        :', h_local)
        print('Seed radius              :', seed_radius)
        print('Initial core nodes       :', len(initial_core))

        interval_time = None
        if SAMPLE_INTERVAL_N is not None and SAMPLE_INTERVAL_N > 0.0:
            interval_time = SAMPLE_INTERVAL_N / CYCLES_PER_TIME
        frames = _select_frames(step.frames, interval_time)
        print('Frames used              : %d / %d' % (len(frames), len(step.frames)))

        rows = []

        for k, frame in enumerate(frames):
            phi_map = _get_phase_field(frame, inst)
            if not phi_map:
                continue

            core_now = set(label for label, phi in phi_map.items()
                           if phi >= TIP_THRESHOLD)
            connect_now = set(label for label, phi in phi_map.items()
                              if phi >= CONNECT_THRESHOLD)

            # Remove the old pre-crack from the traversal graph, except for a
            # small zone around the original tip. This prevents the old crack
            # length from being counted as new crack extension.
            allowed = (connect_now - initial_connect) | seed_band
            allowed.add(root)

            dist = _dijkstra(root, allowed, graph)

            # Only nodes that were not already part of the initial crack can
            # become the new crack tip.
            new_core = core_now - initial_core
            reachable_candidates = [label for label in new_core if label in dist]

            if reachable_candidates:
                tip_label = max(reachable_candidates, key=lambda lab: dist[lab])
                a_geo = float(dist[tip_label])
                tip_xyz = coords[tip_label]
                tip_euclid = _distance(tip_xyz, INITIAL_CRACK_TIP)
            else:
                tip_label = root
                a_geo = 0.0
                tip_xyz = coords[root]
                tip_euclid = 0.0

            N = float(frame.frameValue) * CYCLES_PER_TIME
            x = float(tip_xyz[0]) if len(tip_xyz) > 0 else 0.0
            y = float(tip_xyz[1]) if len(tip_xyz) > 1 else 0.0
            z = float(tip_xyz[2]) if len(tip_xyz) > 2 else 0.0

            rows.append((N, a_geo, x, y, z, tip_euclid,
                         len(new_core), len(reachable_candidates)))

            if k % 25 == 0 or k == len(frames) - 1:
                print('  N=%g, a_geo=%g, tip=(%g, %g)' % (N, a_geo, x, y))

    finally:
        odb.close()

    _ensure_result_dir()
    base = os.path.splitext(os.path.basename(odb_path))[0]
    out_csv = os.path.join(RESULT_DIR, base + '_crack_curve.csv')

    mode = 'wb' if sys.version_info[0] < 3 else 'w'
    with open(out_csv, mode) as f:
        w = csv.writer(f)
        w.writerow([
            'N_cycles',
            'crack_extension_geodesic_mm',
            'tip_x_mm',
            'tip_y_mm',
            'tip_z_mm',
            'tip_straight_distance_mm',
            'n_new_core_nodes',
            'n_reachable_new_core_nodes'
        ])
        for row in rows:
            w.writerow(row)

    print('Wrote:', out_csv)
    return out_csv


def main():
    if len(sys.argv) < 2:
        raise SystemExit(
            'Usage:\n'
            '  abaqus python extract_crack_curve_general.py jobname.odb\n'
            '  abaqus python extract_crack_curve_general.py --all'
        )

    if sys.argv[1].lower() == '--all':
        if not os.path.isdir(ODB_DIR):
            raise SystemExit('ODB_DIR does not exist: %s' % ODB_DIR)
        odb_paths = []
        for name in sorted(os.listdir(ODB_DIR)):
            if name.lower().endswith('.odb'):
                odb_paths.append(os.path.join(ODB_DIR, name))
        if not odb_paths:
            raise SystemExit('No ODB files found in: %s' % ODB_DIR)
    else:
        odb_paths = [_resolve_odb_path(arg) for arg in sys.argv[1:]]

    print('ODB directory   :', ODB_DIR)
    print('Result directory:', RESULT_DIR)
    print('Number of ODBs  :', len(odb_paths))

    failed = []
    for path in odb_paths:
        try:
            extract_one_odb(path)
        except Exception as exc:
            print('FAILED:', path)
            print('Reason:', exc)
            failed.append(path)

    if failed:
        print('\nFailed ODBs:')
        for path in failed:
            print('  ', path)
        raise SystemExit(1)


if __name__ == '__main__':
    main()
