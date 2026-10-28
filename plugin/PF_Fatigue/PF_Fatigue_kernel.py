# -*- coding: utf-8 -*-
# PF-Fatigue v1.0 kernel: model-side Abaqus/CAE operations.
#
# Two simulation strategies, each with the subroutine released with the paper:
#   Cycle-by-cycle : PF_Fatigue_Exact_EnergySplit.f  (triangular waveform, fixed IPC)
#   Envelope (dN)  : PF_Fatigue_GC_Envelope_Final.f  (peak load held, step time = N, dt = dN)

from abaqus import mdb
from abaqusConstants import *
# Attaches the step/output/load constructors to Model objects also in noGUI scripts
# (verified on Abaqus 2023: without it Model has no FieldOutputRequest).
from caeModules import *
import regionToolset
import os

import PF_Fatigue_utils as U

_PLUGIN_DIR = os.path.dirname(os.path.abspath(__file__))
_MODES = {'Cycle-by-cycle': 'cycle', 'Envelope (dN)': 'envelope'}
_MESH_SAMPLE_LIMIT = 60000


def _banner(text):
    print('')
    print('=' * 72)
    print('PF-Fatigue v%s | %s' % (U.VERSION, text))
    print('=' * 72)


def _is_numeric(value):
    try:
        float(value)
        return True
    except Exception:
        return False


def _mode_key(mode):
    if mode not in _MODES:
        raise ValueError('Unknown simulation strategy: %s' % mode)
    return _MODES[mode]


def _collect_section_users(model, material_name):
    users = []
    for sec_name in model.sections.keys():
        try:
            if model.sections[sec_name].material == material_name:
                users.append(sec_name)
        except Exception:
            pass
    return users


def _part_region(part):
    """Whole part: cells (3D), faces (2D) or, for an orphan mesh, all elements."""
    if len(part.cells) > 0:
        return regionToolset.Region(cells=part.cells)
    if len(part.faces) > 0:
        return regionToolset.Region(faces=part.faces)
    return regionToolset.Region(elements=part.elements)


def _auto_section_check(model):
    """Errors/notes for the automatic section: every part must still be without a section."""
    errors, notes = [], []
    parts = sorted(model.parts.keys())
    if not parts:
        return ['The model has no part.'], notes
    # parts that already carry the automatic section (an earlier Build) are fine
    assigned = [p for p in parts
                if any(sa.sectionName != U.AUTO_SECTION_NAME for sa in model.parts[p].sectionAssignments)]
    if assigned:
        errors.append('Part(s) %s already have a section assignment. Choose the material of that section as '
                      'Base material instead of "%s".' % (', '.join(assigned), U.AUTO_SECTION))
    else:
        notes.append('Solid section %s (PF-Fatigue material; thickness 1 for 2D) will be created and assigned to '
                     'part(s): %s.' % (U.AUTO_SECTION_NAME, ', '.join(parts)))
    return errors, notes


def _ensure_auto_section(model, pfMaterialName):
    if U.AUTO_SECTION_NAME not in model.sections:
        model.HomogeneousSolidSection(name=U.AUTO_SECTION_NAME, material=pfMaterialName, thickness=1.0)
    else:
        model.sections[U.AUTO_SECTION_NAME].setValues(material=pfMaterialName)
    done = []
    for p in sorted(model.parts.keys()):
        part = model.parts[p]
        if len(part.sectionAssignments) == 0:
            part.SectionAssignment(region=_part_region(part), sectionName=U.AUTO_SECTION_NAME)
            done.append(p)
    return done


def _instances(model):
    try:
        return [model.rootAssembly.instances[k] for k in model.rootAssembly.instances.keys()]
    except Exception:
        return []


def _element_types(model):
    types = set()
    for inst in _instances(model):
        try:
            for elem in inst.elements:
                types.add(str(elem.type))
        except Exception:
            pass
    return sorted(types)


def _element_check(model):
    types = _element_types(model)
    if not types:
        return [], [], ['Mesh element types could not be checked (model may be unmeshed or use an independent instance without mesh).']
    bad = [t for t in types if t not in U.SUPPORTED_ELEMENTS]
    if bad:
        return ['PF-Fatigue v%s supports %s. Unsupported element type(s): %s'
                % (U.VERSION, ', '.join(U.SUPPORTED_ELEMENTS), ', '.join(bad))], [], []
    return [], ['Element types: %s' % ', '.join(types)], []


def _mesh_size(model, ell):
    """Shortest-edge statistics of the mesh (sampled for very large meshes)."""
    def coords_of(inst):
        nodes = inst.nodes
        elems = inst.elements
        n = len(elems)
        stride = max(1, n // _MESH_SAMPLE_LIMIT)
        for i in range(0, n, stride):
            conn = elems[i].connectivity
            # corner nodes only (first 3/4/8 of the connectivity)
            yield [nodes[c].coordinates for c in conn[:8]]
    elements = []
    sampled = False
    for inst in _instances(model):
        try:
            if len(inst.elements) > _MESH_SAMPLE_LIMIT:
                sampled = True
            elements.extend(coords_of(inst))
        except Exception:
            pass
    h_min, h_med, n = U.element_size_stats(elements)
    notes, warns = U.mesh_advice(ell, h_min, h_med)
    if sampled and h_min is not None:
        notes.append('Mesh statistics from a uniform sample of %d elements.' % n)
    return notes, warns


def _get_load_object(model, source_text):
    kind, name = U.parse_load_source(source_text)
    if kind == 'BC':
        if name not in model.boundaryConditions:
            raise ValueError('Boundary condition not found: %s' % name)
        return kind, name, model.boundaryConditions[name]
    if name not in model.loads:
        raise ValueError('Load not found: %s' % name)
    return kind, name, model.loads[name]


def _states(model, kind, name, stepName):
    """Objects carrying the values of a Load/BC: the object itself, then its step states.

    In Abaqus/CAE the magnitudes of a Load/BC created in an analysis step are stored in
    model.steps[...].loadStates / boundaryConditionStates, not on the object itself.
    The step to be configured is checked first, then all other steps.
    """
    if kind == 'BC':
        obj = model.boundaryConditions[name]
    else:
        obj = model.loads[name]
    yield obj
    steps = [s for s in model.steps.keys() if s != 'Initial']
    if stepName in steps:
        steps.remove(stepName)
        steps.insert(0, stepName)
    for s in steps:
        try:
            states = model.steps[s].boundaryConditionStates if kind == 'BC' else model.steps[s].loadStates
            if name in states.keys():
                yield states[name]
        except Exception:
            pass


def _type_name(obj):
    try:
        return type(obj).__name__
    except Exception:
        return ''


def _creation_step(model, kind, name):
    """First step (in analysis order, Initial included) in which the Load/BC is active.

    Every step lists a state for every Load/BC; steps before its creation carry NOT_YET_ACTIVE
    (verified on Abaqus 2023)."""
    for s in model.steps.keys():
        try:
            states = model.steps[s].boundaryConditionStates if kind == 'BC' else model.steps[s].loadStates
            if name in states.keys() and 'NOT_YET_ACTIVE' not in str(states[name].status):
                return s
        except Exception:
            pass
    return None


def _driver_candidates(model, stepName):
    """Loads/BCs with exactly one non-zero driving component, as 'BC | name (u2 = 0.001)'."""
    out = []
    for kind, repo in (('BC', model.boundaryConditions), ('LOAD', model.loads)):
        for n in sorted(repo.keys()):
            try:
                comp, val, fam = _detect_driver_any(model, kind, n, stepName)
                out.append('%s | %s (%s = %g)' % (kind, n, comp, val))
            except Exception:
                pass
    return out


def _resolve_load_source(model, loadSource, stepName):
    """Empty or invalid selection: use the only Load/BC with one non-zero component, if there is exactly one.
    Returns (load_source_text, note or None)."""
    try:
        kind, name = U.parse_load_source(loadSource)
        repo = model.boundaryConditions if kind == 'BC' else model.loads
        if name in repo.keys():
            return loadSource, None
    except Exception:
        pass
    cands = _driver_candidates(model, stepName)
    if len(cands) == 1:
        chosen = cands[0].split(' (')[0]
        return chosen, ('No valid Load/BC was received from the dialog ("%s"); using the only candidate %s.'
                        % (loadSource, cands[0]))
    return loadSource, None


def _detect_driver_any(model, kind, name, stepName):
    if kind == 'BC' and 'Temperature' in _type_name(model.boundaryConditions[name]):
        raise ValueError('%s is a temperature (phase-field) boundary condition and cannot drive the load. '
                         'Use the Initial crack option for phi = 1 instead.' % name)
    last = None
    for obj in _states(model, kind, name, stepName):
        try:
            return _detect_driver(obj)
        except ValueError as exc:
            last = exc
            if 'several non-zero' in str(exc):
                raise
    raise last


def _detect_driver(obj):
    """Return (component, value, family); exactly one non-zero driving component is required."""
    families = (
        ('displacement', ('u1', 'u2', 'u3', 'ur1', 'ur2', 'ur3')),
        ('force', ('cf1', 'cf2', 'cf3')),
        ('moment', ('cm1', 'cm2', 'cm3')),
    )
    for family, names in families:
        present = []
        for name in names:
            if hasattr(obj, name):
                try:
                    value = getattr(obj, name)
                except Exception:
                    continue
                if _is_numeric(value) and abs(float(value)) > 1.0e-30:
                    present.append((name, float(value)))
        if len(present) == 1:
            return present[0][0], present[0][1], family
        if len(present) > 1:
            raise ValueError('The selected object has several non-zero %s components; '
                             'PF-Fatigue requires one driving component.' % family)
    if hasattr(obj, 'magnitude'):
        try:
            value = getattr(obj, 'magnitude')
            if _is_numeric(value):
                return 'magnitude', float(value), 'scalar'
        except Exception:
            pass
    raise ValueError('Could not detect a non-zero driving component. Supported: displacement/rotation BC, '
                     'concentrated force, moment and scalar (pressure-like) loads.')


def preflight(modelName, baseMaterial, pfMaterialName, E, nu, Gc, ell, split, mode,
              alphaT, kres, etaT, pExp, kappa,
              stepMode, stepName, previousStep, loadSource, maxValue, loadRatio,
              cycles, ipc, deltaN, outputFreq, useInitialCrack, crackSetName,
              jobName, numCpus, memory, workDir):
    """Return (errors, notes, warnings). Nothing in the model is modified."""
    errors, notes, warns = [], [], []
    if modelName not in mdb.models:
        return ['Model not found: %s' % modelName], notes, warns
    model = mdb.models[modelName]
    try:
        key = _mode_key(mode)
    except ValueError as exc:
        return [str(exc)], notes, warns
    notes.append('Simulation strategy: %s' % mode)

    auto = (baseMaterial == U.AUTO_SECTION)
    if auto:
        e_auto, n_auto = _auto_section_check(model)
        errors.extend(e_auto); notes.extend(n_auto)
    elif baseMaterial not in model.materials and pfMaterialName not in model.materials:
        errors.append('Base material "%s" does not exist in the model. Choose the material of your section, '
                      'or "%s" to let PF-Fatigue create and assign the section.' % (baseMaterial, U.AUTO_SECTION))
    if not pfMaterialName:
        errors.append('PF-Fatigue material name is empty.')
    errors.extend(U.validate_material(E, nu, Gc, ell))
    try:
        U.split_code(split)
    except ValueError as exc:
        errors.append(str(exc))

    if key == 'cycle':
        errors.extend(U.validate_cycle_inputs(maxValue, loadRatio, cycles, ipc, alphaT, kres))
        if not errors:
            a0 = U.alpha_t0(Gc, ell)
            notes.append('alpha_T = %g (Gc/(12 l0) = %g).' % (float(alphaT), a0))
    else:
        errors.extend(U.validate_envelope_inputs(maxValue, loadRatio, cycles, deltaN, etaT, pExp, kappa))
        if not errors:
            dns = U.scaled_increment(deltaN, loadRatio, etaT, pExp, kappa)
            notes.append('Scaled cycle increment dN* = dN beta^(2 kappa p)/eta_T = %.4g.' % dns)
            if dns > U.DN_STAR_WARN:
                warns.append('dN* = %.3g exceeds %.1f; the paper uses dN* = 0.625. Check the result with a '
                             'smaller dN.' % (dns, U.DN_STAR_WARN))
            notes.append('Peak load held constant; R = %g enters the material constants only.' % float(loadRatio))
    try:
        if int(outputFreq) < 1:
            errors.append('Field output frequency must be >= 1 increment.')
    except Exception:
        errors.append('Field output frequency must be an integer.')

    if stepMode == 'Create new':
        if previousStep not in model.steps:
            errors.append('Previous step not found: %s' % previousStep)
        if stepName == 'Initial':
            errors.append('Fatigue step cannot be named Initial.')
        if stepName in model.steps:
            errors.append('A step named %s already exists. Choose Step mode = Use existing, or another name.'
                          % stepName)
    elif stepMode == 'Use existing':
        if stepName not in model.steps or stepName == 'Initial':
            errors.append('Existing fatigue step not found: %s. Select an analysis step (not Initial).' % stepName)
        else:
            tn = _type_name(model.steps[stepName])
            if tn.endswith('Step') and 'CoupledTempDisplacement' not in tn:
                errors.append('Step %s is a %s. PF-Fatigue needs a Coupled temp-displacement step: replace it '
                              '(Step Manager -> Replace) or use Step mode = Create new.' % (stepName, tn))
    else:
        errors.append('Unknown step mode: %s' % stepMode)

    loadSource, ls_note = _resolve_load_source(model, loadSource, stepName)
    if ls_note:
        warns.append(ls_note)
    try:
        kind, source_name, source_obj = _get_load_object(model, loadSource)
    except Exception as exc:
        kind = None
        cands = _driver_candidates(model, stepName)
        errors.append(str(exc) + (' Candidates in this model: %s.' % '; '.join(cands) if cands else
                                  ' No Load/BC in this model has exactly one non-zero component.'))
    if kind is not None:
        try:
            comp, old_value, family = _detect_driver_any(model, kind, source_name, stepName)
            notes.append('Load source: %s %s | detected %s = %g' % (kind, source_name, comp, old_value))
        except Exception as exc:
            msg = str(exc)
            if 'Could not detect' in msg:
                msg = ('%s %s has no non-zero component in any step (is it a support?). ' % (kind, source_name)
                       + 'PF-Fatigue needs the Load/BC that drives the specimen, with one non-zero component.')
            cands = [c for c in _driver_candidates(model, stepName) if not c.startswith('%s | %s ' % (kind, source_name))]
            if cands:
                msg += ' Candidates in this model: %s. Select one under Fatigue Loading -> Load / BC.' % '; '.join(cands)
            else:
                msg += ' No Load/BC in this model has exactly one non-zero component.'
            errors.append(msg)
        if stepMode == 'Create new' and previousStep in model.steps:
            order = list(model.steps.keys())
            created = _creation_step(model, kind, source_name)
            if created is not None and order.index(created) > order.index(previousStep):
                errors.append('%s %s is created in step %s, which comes after Previous step %s, so it does not '
                              'exist in the new fatigue step. Use Step mode = Use existing with step %s, or set '
                              'Previous step = %s.' % (kind, source_name, created, previousStep, created, created))

    if bool(useInitialCrack):
        if not crackSetName:
            errors.append('Initial crack is enabled but no Set is selected.')
        elif crackSetName not in model.rootAssembly.sets:
            errors.append('Initial crack Set not found: %s' % crackSetName)

    sec_users = [] if auto else _collect_section_users(model, baseMaterial)
    pf_users = _collect_section_users(model, pfMaterialName)
    if auto:
        pass
    elif not sec_users and not pf_users:
        errors.append('No section uses material "%s". Assign a solid section with this material in the Property '
                      'module, or choose "%s".' % (baseMaterial, U.AUTO_SECTION))
    else:
        if sec_users:
            notes.append('Sections to redirect to the PF-Fatigue material: %s' % ', '.join(sec_users))
        if pf_users:
            notes.append('Sections already using the PF-Fatigue material: %s' % ', '.join(pf_users))

    e_err, e_notes, e_warn = _element_check(model)
    errors.extend(e_err); notes.extend(e_notes); warns.extend(e_warn)
    if not e_err and not U.validate_material(E, nu, Gc, ell):
        m_notes, m_warns = _mesh_size(model, ell)
        notes.extend(m_notes); warns.extend(m_warns)

    try:
        if int(numCpus) < 1:
            errors.append('Number of CPUs must be >= 1.')
        if not (1 <= int(memory) <= 100):
            errors.append('Memory percentage must be between 1 and 100.')
    except Exception:
        errors.append('CPUs and memory must be integers.')
    if not jobName:
        errors.append('Job name is empty.')
    if not workDir:
        errors.append('Working directory is empty.')

    s_err, s_notes = U.check_subroutine(key, _PLUGIN_DIR)
    errors.extend(s_err); notes.extend(s_notes)
    return errors, notes, warns


def _ensure_pf_material(model, baseMaterial, pfMaterialName, mech):
    if pfMaterialName not in model.materials:
        mat = model.Material(name=pfMaterialName,
                             description='PF-Fatigue v%s | phase-field fatigue (UMAT + UMATHT)' % U.VERSION)
    else:
        mat = model.materials[pfMaterialName]
    try:
        mat.depvar.setValues(n=4)
    except Exception:
        mat.Depvar(n=4)
    # A thermomechanical user material supplies the UMAT constants and the
    # (unused) UMATHT constant. Density is an interface placeholder only.
    try:
        mat.userMaterial.setValues(type=THERMOMECHANICAL, mechanicalConstants=mech, thermalConstants=(1.0,))
    except Exception:
        mat.UserMaterial(type=THERMOMECHANICAL, mechanicalConstants=mech, thermalConstants=(1.0,))
    try:
        mat.density
    except Exception:
        mat.Density(table=((1.0,),))
    redirected = []
    for sec_name in model.sections.keys():
        sec = model.sections[sec_name]
        try:
            if sec.material == baseMaterial and baseMaterial != pfMaterialName:
                sec.setValues(material=pfMaterialName)
                redirected.append(sec_name)
        except Exception:
            pass
    return redirected


def _ensure_step(model, stepMode, stepName, previousStep, timePeriod, dt, maxNumInc):
    # Fixed increments are enforced by initial = minimum = maximum increment, as in the article decks
    # ("dt, period, dt, dt"). timeIncrementationMethod=FIXED alone is NOT honoured for a steady-state
    # coupled step in Abaqus 2023 (the data line is written as "dt, period," and Abaqus grows dt).
    kwargs = dict(response=STEADY_STATE, timePeriod=float(timePeriod), nlgeom=OFF,
                  timeIncrementationMethod=AUTOMATIC, maxNumInc=int(maxNumInc),
                  initialInc=float(dt), minInc=float(dt), maxInc=float(dt), creepIntegration=None,
                  solutionTechnique=SEPARATED, matrixStorage=SYMMETRIC)

    def apply(kw):
        if stepMode == 'Create new' and stepName not in model.steps:
            model.CoupledTempDisplacementStep(name=stepName, previous=previousStep, **kw)
        else:
            model.steps[stepName].setValues(**kw)
    try:
        apply(kwargs)               # creepIntegration=None writes "creep=none"
    except TypeError:
        kwargs.pop('creepIntegration')
        apply(kwargs)


def _ensure_amplitude(model, amplitudeName, data):
    if amplitudeName in model.amplitudes:
        model.amplitudes[amplitudeName].setValues(data=data, smooth=0.0, timeSpan=STEP)
    else:
        model.TabularAmplitude(name=amplitudeName, data=data, smooth=0.0, timeSpan=STEP)


def _apply_amplitude_to_source(model, stepName, loadSource, maxValue, amplitudeName):
    kind, name, obj = _get_load_object(model, loadSource)
    component, old_value, family = _detect_driver_any(model, kind, name, stepName)
    kwargs = {'amplitude': amplitudeName, component: float(maxValue)}
    try:
        obj.setValuesInStep(stepName=stepName, **kwargs)
    except Exception as exc:
        raise ValueError('Failed to apply the amplitude to %s "%s" (component %s): %s'
                         % (kind, name, component, str(exc)))
    return kind, name, component


def _ensure_output(model, stepName, outputFreq):
    name = 'PF-Fatigue-F-Output'
    variables = ('U', 'RF', 'S', 'NT', 'SDV')
    # The repository exists only after a first output request has been defined.
    existing = getattr(model, 'fieldOutputRequests', None)
    if existing is not None and name in existing.keys():
        model.fieldOutputRequests[name].setValuesInStep(stepName=stepName, variables=variables,
                                                        frequency=int(outputFreq))
    else:
        model.FieldOutputRequest(name=name, createStepName=stepName, variables=variables,
                                 frequency=int(outputFreq))
    return name


def _ensure_initial_crack(model, stepName, useInitialCrack, crackSetName, timePeriod):
    """Initial crack: phi = 1 held on the Set by a temperature BC with a constant amplitude.

    An initial temperature (initial condition) is not enough: the step is steady state, so the phase-field
    equation does not remember an initial value. A BC without amplitude would be ramped over the step
    (phi = t/T), so a constant amplitude is attached. Abaqus does not accept a non-zero BC in the Initial
    step; the BC is created in the fatigue step."""
    bc_name, amp_name = 'PF-Fatigue-InitialCrack', 'PF-Fatigue-Const'
    if not bool(useInitialCrack):
        # an initial-crack BC left by an earlier Build of this plugin must not stay active
        if bc_name in model.boundaryConditions.keys():
            del model.boundaryConditions[bc_name]
            return 'Skipped (earlier BC %s deleted)' % bc_name
        return 'Skipped'
    if bc_name in model.predefinedFields.keys():          # v1.0 pre-release used an initial temperature
        del model.predefinedFields[bc_name]
    data = ((0.0, 1.0), (float(timePeriod), 1.0))
    if amp_name in model.amplitudes.keys():
        model.amplitudes[amp_name].setValues(data=data, smooth=0.0, timeSpan=STEP)
    else:
        model.TabularAmplitude(name=amp_name, data=data, smooth=0.0, timeSpan=STEP)
    region = model.rootAssembly.sets[crackSetName]
    if bc_name in model.boundaryConditions.keys():
        del model.boundaryConditions[bc_name]
    model.TemperatureBC(name=bc_name, createStepName=stepName, region=region, magnitude=1.0,
                        amplitude=amp_name, distributionType=UNIFORM)
    return '%s (phi = 1 held by temperature BC %s)' % (crackSetName, bc_name)


def _ensure_job(modelName, jobName, subroutine, numCpus, memory, workDir, submitJob):
    if not os.path.isdir(workDir):
        os.makedirs(workDir)
    os.chdir(workDir)
    if jobName in mdb.jobs:
        del mdb.jobs[jobName]
    # Abaqus/Standard requires the number of domains to be a multiple of the number of CPUs.
    job = mdb.Job(name=jobName, model=modelName, type=ANALYSIS, userSubroutine=subroutine,
                  numCpus=int(numCpus), numDomains=int(numCpus), memory=int(memory), memoryUnits=PERCENTAGE)
    if bool(submitJob):
        job.submit(consistencyChecking=OFF)


def configure_and_run(modelName, baseMaterial, pfMaterialName, E, nu, Gc, ell, split, mode,
                      alphaT, kres, etaT, pExp, kappa,
                      stepMode, stepName, previousStep, loadSource, maxValue, loadRatio,
                      cycles, ipc, deltaN, amplitudeName, outputFreq,
                      useInitialCrack, crackSetName,
                      jobName, numCpus, memory, workDir, submitJob, checkOnly=False):
    _banner('Check Model' if checkOnly else 'Build / Update')
    print('[INFO] Received: model=%r, baseMaterial=%r, stepMode=%r, step=%r, loadSource=%r'
          % (modelName, baseMaterial, stepMode, stepName, loadSource))
    if modelName in mdb.models:
        loadSource, ls_note = _resolve_load_source(mdb.models[modelName], loadSource, stepName)
        if ls_note:
            print('[WARN] ' + ls_note)
    errors, notes, warns = preflight(
        modelName, baseMaterial, pfMaterialName, E, nu, Gc, ell, split, mode,
        alphaT, kres, etaT, pExp, kappa, stepMode, stepName, previousStep, loadSource,
        maxValue, loadRatio, cycles, ipc, deltaN, outputFreq, useInitialCrack, crackSetName,
        jobName, numCpus, memory, workDir)
    for item in notes:
        print('[OK]   ' + item)
    for item in warns:
        print('[WARN] ' + item)
    if errors:
        print('')
        print('Preflight FAILED:')
        for item in errors:
            print('  [ERROR] ' + item)
        raise ValueError('PF-Fatigue preflight failed. See the message area.')
    if bool(checkOnly):
        print('')
        print('Preflight PASSED. No model objects were modified.')
        print('=' * 72)
        return

    model = mdb.models[modelName]
    key = _mode_key(mode)
    if key == 'cycle':
        mech = U.cycle_props(E, nu, Gc, ell, alphaT, kres, split)
        amp, period, dt, total = U.build_triangle_amplitude(cycles, loadRatio, ipc)
    else:
        mech = U.envelope_props(E, nu, Gc, ell, etaT, split, pExp, loadRatio, kappa)
        amp, period, dt, total = U.build_envelope_schedule(cycles, deltaN)
    subroutine = U.subroutine_path(key, _PLUGIN_DIR)

    redirected = _ensure_pf_material(model, baseMaterial, pfMaterialName, mech)
    if baseMaterial == U.AUTO_SECTION:
        redirected = ['%s -> part(s) %s' % (U.AUTO_SECTION_NAME, ', '.join(_ensure_auto_section(model, pfMaterialName)))]
    _ensure_step(model, stepMode, stepName, previousStep, period, dt, total + 10)
    _ensure_amplitude(model, amplitudeName, amp)
    kind, src_name, component = _apply_amplitude_to_source(model, stepName, loadSource, maxValue, amplitudeName)
    out_name = _ensure_output(model, stepName, outputFreq)
    crack_status = _ensure_initial_crack(model, stepName, useInitialCrack, crackSetName, period)
    _ensure_job(modelName, jobName, subroutine, numCpus, memory, workDir, submitJob)

    print('')
    print('PF-Fatigue configuration completed.')
    print('  Strategy           : %s' % mode)
    print('  Material           : %s  PROPS = %s' % (pfMaterialName, ', '.join('%g' % v for v in mech)))
    if redirected:
        print('  Redirected sections: %s' % ', '.join(redirected))
    print('  Fatigue step       : %s  (period %g, fixed increment %g, %d increments)' % (stepName, period, dt, total))
    if key == 'cycle':
        print('  Waveform           : triangular from zero, max %g, R = %g, IPC = %d'
              % (float(maxValue), float(loadRatio), int(ipc)))
    else:
        print('  Loading            : peak %g held; R = %g in PROPS; dN = %g; dN* = %.4g'
              % (float(maxValue), float(loadRatio), float(deltaN),
                 U.scaled_increment(deltaN, loadRatio, etaT, pExp, kappa)))
    print('  Load source        : %s | %s | %s' % (kind, src_name, component))
    print('  Field output       : %s every %d increments' % (out_name, int(outputFreq)))
    print('  Initial crack      : %s' % crack_status)
    print('  Job                : %s  (%s)' % (jobName, os.path.basename(subroutine)))
    print('  Working directory  : %s' % os.path.abspath(workDir))
    print('  Submitted          : %s' % bool(submitJob))
    print('=' * 72)
