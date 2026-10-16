# -*- coding: utf-8 -*-
# PF-Fatigue v1.0 GUI. User-facing controls only; PROPS and SDVs stay internal.
# Compact two-column layout: the dialog must fit on a 1080-pixel screen.

from abaqusConstants import *
from abaqusGui import *
from kernelAccess import mdb
import os
import PF_Fatigue_utils as U


def _note(parent, text):
    FXLabel(p=parent, text=text, opts=JUSTIFY_LEFT)


def _group(parent, title, fillY=False):
    opts = FRAME_GROOVE | LAYOUT_FILL_X | (LAYOUT_FILL_Y if fillY else 0)
    return FXGroupBox(parent, title, opts)


def _columns(parent):
    hf = FXHorizontalFrame(p=parent, opts=LAYOUT_FILL_X | LAYOUT_FILL_Y, pl=0, pr=0, pt=0, pb=0)
    left = FXVerticalFrame(p=hf, opts=LAYOUT_FILL_X | LAYOUT_FILL_Y, pl=0, pr=0, pt=0, pb=0)
    right = FXVerticalFrame(p=hf, opts=LAYOUT_FILL_X | LAYOUT_FILL_Y, pl=0, pr=0, pt=0, pb=0)
    return left, right


def _field(parent, label, kw, ncols=12):
    return AFXTextField(p=parent, ncols=ncols, labelText=label, tgt=kw, sel=0)


def _combo(parent, label, kw, items, ncols=20):
    cb = AFXComboBox(p=parent, ncols=ncols, nvis=1, text=label, tgt=kw, sel=0)
    cb.setMaxVisible(15)
    for it in items:
        cb.appendItem(text=it)
    return cb


class PF_FatigueDB(AFXDataDialog):

    ID_CHECK_MODEL = AFXDataDialog.ID_LAST
    ID_WORKDIR = AFXDataDialog.ID_LAST + 1

    def __init__(self, form):
        AFXDataDialog.__init__(self, form, 'PF-Fatigue v1.0',
                               self.APPLY | self.CANCEL, DIALOG_ACTIONS_SEPARATOR)
        self.form = form
        self.dirDb = None
        self.logoIcon = None
        self._state = None
        # Build / Update is the Apply button, so the dialog stays open after building.
        self.getActionButton(self.ID_CLICKED_APPLY).setText('Build / Update')
        self.getActionButton(self.ID_CLICKED_CANCEL).setText('Close')
        self.appendActionButton('Check Model', self, self.ID_CHECK_MODEL)
        FXMAPFUNC(self, SEL_COMMAND, self.ID_CHECK_MODEL, AFXDataDialog.onCmdApply)

        main = FXVerticalFrame(p=self, opts=LAYOUT_FILL_X | LAYOUT_FILL_Y, pl=8, pr=8, pt=6, pb=4)

        # Header: two short lines and a small logo
        header = FXHorizontalFrame(p=main, opts=LAYOUT_FILL_X, pl=0, pr=0, pt=0, pb=2)
        headerText = FXVerticalFrame(p=header, opts=LAYOUT_FILL_X | LAYOUT_CENTER_Y, pl=0, pr=0, pt=0, pb=0)
        title = FXLabel(p=headerText, text='Gradient-consistent phase-field fatigue (UMAT + UMATHT)',
                        opts=JUSTIFY_LEFT)
        title.setFont(getAFXFont(FONT_BOLD))
        _note(headerText, 'Use one consistent unit system (e.g. N, mm, MPa). '
                          'Reopen the dialog after changing the active model.')
        try:
            self.logoIcon = afxCreatePNGIcon(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                          'icons', 'xjtu_logo_44.png'))
            FXLabel(header, '', self.logoIcon, opts=LAYOUT_RIGHT | LAYOUT_CENTER_Y)
        except Exception:
            pass
        tabs = FXTabBook(p=main, tgt=None, sel=0, opts=TABBOOK_NORMAL | LAYOUT_FILL_X | LAYOUT_FILL_Y)

        # Current model inventory
        modelNames = sorted(list(mdb.models.keys()))
        currentModel = ''
        try:
            currentModel = getAFXApp().getAFXMainWindow().getCurrentContext().get('modelName', '')
        except Exception:
            currentModel = ''
        if currentModel not in modelNames and modelNames:
            currentModel = modelNames[0]
        if currentModel:
            form.modelNameKw.setValue(currentModel)
        model = mdb.models[currentModel] if currentModel in mdb.models else None
        form.stepModeKw.setValue('Create new')
        form.stepNameKw.setValue('PF-Fatigue-Step')

        # ------------------------------------------------------------- MODEL
        FXTabItem(p=tabs, text='Model')
        t = FXVerticalFrame(p=tabs, opts=FRAME_RAISED | LAYOUT_FILL_X | LAYOUT_FILL_Y, pl=6, pr=6, pt=4, pb=4)
        left, right = _columns(t)
        va = AFXVerticalAligner(_group(left, 'Model and material', True))
        _combo(va, 'Abaqus model:', form.modelNameKw, modelNames)
        materialNames = sorted(list(model.materials.keys())) if model is not None else []
        _combo(va, 'Base material:', form.baseMaterialKw, [U.AUTO_SECTION] + materialNames, 22)
        form.baseMaterialKw.setValue(U.guess_base_material(model, materialNames) if model is not None
                                     else U.AUTO_SECTION)
        _field(va, 'PF-Fatigue material:', form.pfMaterialNameKw, 20)
        _combo(va, 'Energy split:', form.splitKw, ['None', 'Amor', 'Miehe', 'No-tension'], 12)
        va = AFXVerticalAligner(_group(right, 'Elasticity and fracture', True))
        _field(va, 'Young modulus E:', form.EKw)
        _field(va, 'Poisson ratio nu:', form.nuKw)
        _field(va, 'Fracture energy Gc:', form.GcKw)
        _field(va, 'Length scale l0:', form.ellKw)

        g = _group(t, 'Fatigue law')
        _combo(g, 'Strategy:', form.modeKw, ['Cycle-by-cycle', 'Envelope (dN)'], 16)
        left, right = _columns(g)
        gc = _group(left, 'Cycle-by-cycle (every cycle resolved)', True)
        va = AFXVerticalAligner(gc)
        self.cycleWidgets = [_field(va, 'Threshold alpha_T:', form.alphaTKw),
                             _field(va, 'Residual stiffness k_res:', form.kresKw)]
        _note(gc, 'Common choice: alpha_T = Gc / (12 l0)')
        ge = _group(right, 'Envelope (peak held, dN cycles/incr.)', True)
        va = AFXVerticalAligner(ge)
        self.envelopeWidgets = [_field(va, 'Life scale eta_T:', form.etaTKw),
                                _field(va, 'Exponent p:', form.pExpKw),
                                _field(va, 'Mean-stress kappa:', form.kappaKw)]
        _note(ge, 'p: curve slope; eta_T, kappa: cycle scale')

        # ---------------------------------------------------------- LOADING
        FXTabItem(p=tabs, text='Fatigue Loading')
        t = FXVerticalFrame(p=tabs, opts=FRAME_RAISED | LAYOUT_FILL_X | LAYOUT_FILL_Y, pl=6, pr=6, pt=4, pb=4)
        left, right = _columns(t)

        va = AFXVerticalAligner(_group(left, 'Step'))
        _combo(va, 'Step mode:', form.stepModeKw, ['Create new', 'Use existing'], 18)
        stepNames = list(model.steps.keys()) if model is not None else []
        _combo(va, 'Fatigue step:', form.stepNameKw,
               ['PF-Fatigue-Step'] + [s for s in stepNames if s not in ('Initial', 'PF-Fatigue-Step')], 18)
        self.previousStepCombo = _combo(va, 'Previous step:', form.previousStepKw, stepNames, 18)
        if 'Initial' in stepNames:
            form.previousStepKw.setValue('Initial')
        elif stepNames:
            form.previousStepKw.setValue(stepNames[-1])
        # Default to the step that already holds a driving Load/BC (the usual beginner workflow)
        drv, drvStep = U.guess_driver(model) if model is not None else (None, None)
        if drvStep:
            form.stepModeKw.setValue('Use existing')
            form.stepNameKw.setValue(drvStep)
            form.previousStepKw.setValue(drvStep)

        g = _group(left, 'Load source')
        items = []
        if model is not None:
            items = ['BC | ' + n for n in sorted(model.boundaryConditions.keys())] + \
                    ['LOAD | ' + n for n in sorted(model.loads.keys())]
        cb = _combo(g, 'Load / BC:', form.loadSourceKw, items, 22)
        if items:
            form.loadSourceKw.setValue(drv if drv in items else items[0])
        else:
            cb.appendItem(text='-- No Load / BC found --')
            form.loadSourceKw.setValue('')
        _note(g, 'One non-zero component; PF-Fatigue sets its amplitude.')

        g = _group(left, 'Initial crack (optional)', True)
        FXCheckButton(p=g, text='Set phase field = 1 on an assembly Set', tgt=form.useInitialCrackKw, sel=0)
        setNames = []
        if model is not None:
            try:
                setNames = sorted(list(model.rootAssembly.sets.keys()))
            except Exception:
                setNames = []
        self.crackSetCombo = _combo(g, 'Set:', form.crackSetNameKw, setNames, 22)
        if setNames:
            form.crackSetNameKw.setValue(setNames[0])
        else:
            self.crackSetCombo.appendItem(text='-- No assembly Set found --')
            form.crackSetNameKw.setValue('')

        g = _group(right, 'Cycles', True)
        va = AFXVerticalAligner(g)
        _field(va, 'Maximum (peak) value:', form.maxValueKw)
        _field(va, 'Load ratio R = min/max:', form.loadRatioKw)
        _field(va, 'Number of cycles:', form.cyclesKw)
        self.cycleWidgets.append(_field(va, 'IPC, increments/cycle:', form.ipcKw))
        self.envelopeWidgets.append(_field(va, 'dN, cycles/increment:', form.deltaNKw))
        _field(va, 'Amplitude name:', form.amplitudeNameKw, 16)
        _field(va, 'Output every N incr.:', form.outputFreqKw)
        _note(g, 'Cycle-by-cycle: triangle from zero,\n'
                 '  dt = 1/IPC (IPC a multiple of 4).\n'
                 'Envelope: peak held, dt = dN;\n'
                 '  Check Model reports dN*.')

        # ------------------------------------------------------------ CHECK
        FXTabItem(p=tabs, text='Check')
        t = FXVerticalFrame(p=tabs, opts=FRAME_RAISED | LAYOUT_FILL_X | LAYOUT_FILL_Y, pl=6, pr=6, pt=4, pb=4)
        g = _group(t, 'Check Model (changes nothing)', True)
        _note(g, '- model, base material and section assignment\n'
                 '- element types (CPE4T, CPS4T, CPE3T, C3D8T)\n'
                 '- element size relative to l0 (h <= l0/2 recommended)\n'
                 '- Load / BC with one driving component\n'
                 '- cycle definition; envelope: scaled increment dN*\n'
                 '- optional initial-crack Set, job settings\n'
                 '- bundled subroutine of the selected strategy (SHA256)\n\n'
                 'Messages appear in the message area: [OK], [WARN], [ERROR].')

        # -------------------------------------------------------------- RUN
        FXTabItem(p=tabs, text='Run')
        t = FXVerticalFrame(p=tabs, opts=FRAME_RAISED | LAYOUT_FILL_X | LAYOUT_FILL_Y, pl=6, pr=6, pt=4, pb=4)
        g = _group(t, 'Job', True)
        va = AFXVerticalAligner(g)
        _field(va, 'Job name:', form.jobNameKw, 24)
        _field(va, 'CPUs:', form.numCpusKw, 8)
        _field(va, 'Memory (%):', form.memoryKw, 8)
        row = FXHorizontalFrame(p=g, opts=LAYOUT_FILL_X, pl=0, pr=0, pt=0, pb=0)
        _field(row, 'Working directory:', form.workDirKw, 36)
        FXButton(p=row, text='Browse...', tgt=self, sel=self.ID_WORKDIR)
        FXMAPFUNC(self, SEL_COMMAND, self.ID_WORKDIR, PF_FatigueDB.onCmdWorkDir)
        FXCheckButton(p=g, text='Submit job immediately after Build / Update', tgt=form.submitJobKw, sel=0)
        _note(g, 'The user subroutine is chosen from the strategy (bundled, unmodified).\n'
                 'First use: leave Submit unchecked, inspect the created objects, then submit.')

    def processUpdates(self):
        # Grey out the inputs that the selected strategy / options do not use.
        try:
            state = (self.form.modeKw.getValue(), self.form.stepModeKw.getValue(),
                     bool(self.form.useInitialCrackKw.getValue()))
            if state != self._state:
                self._state = state
                envelope = state[0].startswith('Envelope')
                for w in self.cycleWidgets:
                    w.disable() if envelope else w.enable()
                for w in self.envelopeWidgets:
                    w.enable() if envelope else w.disable()
                self.previousStepCombo.enable() if state[1] == 'Create new' else self.previousStepCombo.disable()
                self.crackSetCombo.enable() if state[2] else self.crackSetCombo.disable()
        except Exception:
            pass
        base = getattr(AFXDataDialog, 'processUpdates', None)
        return base(self) if base is not None else None

    def onCmdWorkDir(self, sender, sel, ptr):
        if not self.dirDb:
            self.dirDb = AFXFileSelectorDialog(self, 'Select PF-Fatigue Working Directory',
                                               self.form.workDirKw, None, AFXSELECTFILE_DIRECTORY)
            self.dirDb.create()
        self.dirDb.showModal()
        return 1
