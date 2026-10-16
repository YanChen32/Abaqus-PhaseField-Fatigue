# -*- coding: utf-8 -*-
# PF-Fatigue v1.0 - Abaqus/CAE GUI registration and command keywords.

from abaqusGui import *
from abaqusConstants import ALL
import os


class PF_Fatigue_plugin(AFXForm):

    def __init__(self, owner):
        AFXForm.__init__(self, owner)
        self.cmd = AFXGuiCommand(mode=self, method='configure_and_run',
                                 objectName='PF_Fatigue_kernel', registerQuery=False)
        c = self.cmd
        # Model and material
        self.modelNameKw = AFXStringKeyword(c, 'modelName', True, 'Model-1')
        self.baseMaterialKw = AFXStringKeyword(c, 'baseMaterial', True, 'Material-1')
        self.pfMaterialNameKw = AFXStringKeyword(c, 'pfMaterialName', True, 'PF-Fatigue-Material')
        self.EKw = AFXFloatKeyword(c, 'E', True, 210000.0)
        self.nuKw = AFXFloatKeyword(c, 'nu', True, 0.30)
        self.GcKw = AFXFloatKeyword(c, 'Gc', True, 2.70)
        self.ellKw = AFXFloatKeyword(c, 'ell', True, 0.004)
        self.splitKw = AFXStringKeyword(c, 'split', True, 'None')
        self.modeKw = AFXStringKeyword(c, 'mode', True, 'Cycle-by-cycle')
        # Cycle-by-cycle law
        self.alphaTKw = AFXFloatKeyword(c, 'alphaT', True, 56.25)
        self.kresKw = AFXFloatKeyword(c, 'kres', True, 1.0e-7)
        # Envelope law
        self.etaTKw = AFXFloatKeyword(c, 'etaT', True, 1.0)
        self.pExpKw = AFXFloatKeyword(c, 'pExp', True, 1.2)
        self.kappaKw = AFXFloatKeyword(c, 'kappa', True, 1.0)
        # Step and loading
        self.stepModeKw = AFXStringKeyword(c, 'stepMode', True, 'Create new')
        self.stepNameKw = AFXStringKeyword(c, 'stepName', True, 'PF-Fatigue-Step')
        self.previousStepKw = AFXStringKeyword(c, 'previousStep', True, 'Initial')
        self.loadSourceKw = AFXStringKeyword(c, 'loadSource', True, '')
        self.maxValueKw = AFXFloatKeyword(c, 'maxValue', True, 0.002)
        self.loadRatioKw = AFXFloatKeyword(c, 'loadRatio', True, -1.0)
        self.cyclesKw = AFXFloatKeyword(c, 'cycles', True, 350.0)
        self.ipcKw = AFXIntKeyword(c, 'ipc', True, 128)
        self.deltaNKw = AFXFloatKeyword(c, 'deltaN', True, 10.0)
        self.amplitudeNameKw = AFXStringKeyword(c, 'amplitudeName', True, 'PF-Fatigue-Amp')
        self.outputFreqKw = AFXIntKeyword(c, 'outputFreq', True, 10)
        # Optional initial phase-field crack
        self.useInitialCrackKw = AFXBoolKeyword(c, 'useInitialCrack', AFXBoolKeyword.TRUE_FALSE, True, False)
        self.crackSetNameKw = AFXStringKeyword(c, 'crackSetName', True, '')
        # Job
        self.jobNameKw = AFXStringKeyword(c, 'jobName', True, 'PF-Fatigue-Job')
        self.numCpusKw = AFXIntKeyword(c, 'numCpus', True, 4)
        self.memoryKw = AFXIntKeyword(c, 'memory', True, 90)
        try:
            defaultWorkDir = getAFXApp().getAFXMainWindow().getWorkDirectory()
        except Exception:
            defaultWorkDir = os.getcwd()
        self.workDirKw = AFXStringKeyword(c, 'workDir', True, defaultWorkDir)
        self.submitJobKw = AFXBoolKeyword(c, 'submitJob', AFXBoolKeyword.TRUE_FALSE, True, False)
        self.checkOnlyKw = AFXBoolKeyword(c, 'checkOnly', AFXBoolKeyword.TRUE_FALSE, True, False)

    def getFirstDialog(self):
        import PF_FatigueDB
        return PF_FatigueDB.PF_FatigueDB(self)

    def doCustomChecks(self):
        checkOnly = False
        try:
            checkOnly = (self.getPressedButtonId() == self.getCurrentDialog().ID_CHECK_MODEL)
        except Exception:
            checkOnly = False
        self.checkOnlyKw.setValue(checkOnly)
        try:
            wd = self.workDirKw.getValue()
            if wd and os.path.isdir(wd):
                getAFXApp().getAFXMainWindow().setWorkDirectory(wd)
        except Exception:
            pass
        return True

    def okToCancel(self):
        return False


thisDir = os.path.dirname(os.path.abspath(__file__))
menuIcon = None
try:
    menuIcon = afxCreatePNGIcon(os.path.join(thisDir, 'icons', 'pf_fatigue_16.png'))
except Exception:
    menuIcon = None

toolset = getAFXApp().getAFXMainWindow().getPluginToolset()
toolset.registerGuiMenuButton(
    buttonText='PF-Fatigue',
    object=PF_Fatigue_plugin(toolset),
    messageId=AFXMode.ID_ACTIVATE,
    icon=menuIcon,
    kernelInitString='import PF_Fatigue_kernel',
    applicableModules=ALL,
    version='1.0',
    author='PF-Fatigue',
    description='Gradient-consistent phase-field fatigue for Abaqus (UMAT + UMATHT)',
    helpUrl='')
