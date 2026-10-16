# PF-Fatigue-GC

Gradient-consistent phase-field fatigue fracture in Abaqus/Standard with a UMAT + UMATHT pair.

The package belongs to the article *A Gradient-Consistent Abaqus UMAT/UMATHT Framework for Phase-Field
Fatigue Fracture* (K. Kang, Y. Chen, Y. Liu; International Journal of Fatigue, submitted). It contains the
user subroutines of the article, the input files of the two SENT examples, scripts for post-processing and
parameter identification, and the Abaqus/CAE plug-in PF-Fatigue v1.0.

The full description is `PF-Fatigue-GC_Documentation.pdf`: implementation,
model set-up in Abaqus/CAE, a verified example, parameter identification and theory. This README is a
summary.

## 1. Principle

The phase field φ is the temperature degree of freedom of the Abaqus coupled temperature–displacement
elements. The fatigue degradation function *f* multiplies the fracture resistance. In the
gradient-consistent (GC) form the phase-field equation is
div(*f* ∇φ) + 2(1−φ)*H*/(*G*c *l*₀) − *f*φ/*l*₀² = 0. UMATHT returns the heat flux q = −*f* ∇φ, so the
element assembly contains the term (∇ ln *f*)·∇φ without a user element. The local-equivalent (LE) form
replaces this by a locally scaled toughness and is included for comparison.

## 2. Contents

| Folder | Content |
|---|---|
| `src/` | Fortran user subroutines (UMAT and UMATHT in one file each) |
| `examples/` | input files of the Mode-I and shear-loaded SENT examples |
| `postprocessing/abaqus_python/` | ODB extraction scripts (run with `abaqus python`) |
| `postprocessing/identification/` | S–N and Paris parameter identification (Python 3, numpy, scipy) and the experimental data of the article |
| `plugin/` | Abaqus/CAE plug-in PF-Fatigue v1.0 and its User Guide (`PF-Fatigue_v1.0_User_Guide.pdf`) |
| `PF-Fatigue-GC_Documentation.pdf` | documentation |
| `LICENSE`, `CITATION.cff` | licence, citation |

## 3. User subroutines

| File | Fatigue law | Use |
|---|---|---|
| `PF_Fatigue_Exact_EnergySplit.f` | cycle by cycle, GC | resolved load cycles (SENT examples) |
| `PF_Fatigue_LE_EnergySplit.f` | cycle by cycle, LE | comparison; differs from the GC file only in the source term and the flux |
| `PF_Fatigue_GC_Envelope_Final.f` | cycle domain (peak-load envelope), GC | long lives, S–N, Paris and structural examples; each increment represents ΔN cycles |

### Material constants (`*User Material, type=MECHANICAL`)

| PROPS | Cycle by cycle (7) | Cycle domain (9) |
|---|---|---|
| 1 | E | E |
| 2 | ν | ν |
| 3 | *G*c | *G*c |
| 4 | *l*₀ | *l*₀ |
| 5 | α_T (fatigue threshold) | η_T (dimensionless threshold) |
| 6 | *k*res (residual stiffness) | ISPLIT |
| 7 | ISPLIT | *p* (exponent) |
| 8 | – | *R* (load ratio) |
| 9 | – | κ (mean-stress parameter) |

ISPLIT selects the energy split of the crack driving force: 0 none, 1 Amor (volumetric–deviatoric),
2 Miehe (spectral), 3 no-tension. The whole stress is degraded (hybrid formulation). The cycle-domain file
uses *k*res = 10⁻⁶ internally. A thermal user material with one constant (value 1) and `*Depvar` = 4 are
also needed. The two kinds of files use different constants and state variables.

### State variables

| SDV | Cycle by cycle | Cycle domain |
|---|---|---|
| 1 | accumulated fatigue variable ᾱ | accumulated history Ā |
| 2 | history field *H* | history field *H* |
| 3 | α of the previous increment | equivalent driving force *A*eq |
| 4 | degradation *f* | degradation *f* |

### Analysis settings

* `*Coupled Temperature-displacement, steady state` with `*Solution Technique, type=SEPARATED`
  (staggered update). The temperature NT11 is the phase field: 0 intact, 1 broken.
* Cycle by cycle: triangular load amplitude; fixed increment 1/IPC of a cycle.
* Cycle domain: the peak load is held; the step time is the cycle number and the fixed increment is ΔN.
  *R* enters only through the material constants. The increment is not adapted during the analysis.
* Tested with Abaqus/Standard 2023 and Intel Fortran 19.1 (Visual Studio 2019) on Windows.

## 4. Quick start

```bat
cd examples\SENT_I_tension
abaqus job=SENT_I_GC input=SENT_I_IPC128.inp user=..\..\src\PF_Fatigue_Exact_EnergySplit.f cpus=8 interactive
abaqus job=SENT_I_LE input=SENT_I_IPC128.inp user=..\..\src\PF_Fatigue_LE_EnergySplit.f cpus=8 interactive
abaqus python ..\..\postprocessing\abaqus_python\extract_crack_curve_SENT.py SENT_I_GC.odb
```

Both examples use N, mm and MPa; run each in its own working directory.

* `examples/SENT_I_tension/SENT_I_IPC128.inp` (article Figs. 2 and 3): IPC = 128, 350 cycles, ISPLIT = 0. For
  another IPC, set the fixed increment (line after `*Coupled Temperature-displacement`) to 1/IPC and the output
  frequency to IPC. Crack extension: `extract_crack_curve_SENT.py`.
* `examples/SENT_II_shear/SENT_II_IPC16.inp` (article Figs. 4 to 6): IPC = 16, 2000 cycles, ISPLIT = 2. The file
  declares `*Depvar` = 7, as in the runs of the article; the subroutines use SDV1–SDV4 only. Crack extension along
  the curved crack: `extract_crack_curve_general.py`; paths at matched crack extension: `export_sent2_path_raw.py`;
  initiation angle: `extract_SENT2_initiation_angles.py`.

Each run can be repeated with `src/PF_Fatigue_LE_EnergySplit.f` for the LE comparison. The input files of the
other examples of the article are not part of this package. The scripts in `postprocessing/abaqus_python/` are
run with `abaqus python <script> <odb>`, open the ODB read-only and write CSV files in the working directory;
settings such as the initial crack-tip position are at the top of each script.

## 5. Numerical settings

* Mesh: element size h ≤ *l*₀/2 along the crack path.
* Cycle by cycle: the result depends on IPC; compare formulations at the same IPC.
* Cycle domain: the resolution is set by the scaled increment ΔN* = ΔN·S^(2κp)/η_T with S = (1−R)/2, not
  by ΔN. A fixed ΔN gives a different resolution at each *R*. The crack-growth calculations of the
  article use ΔN* = 0.625. The plug-in reports ΔN* and warns above 1.
* Initial transient: growth rates in the first ≈20 *l*₀ from a sharp notch are not in the Paris regime
  and are excluded from the identification.

## 6. Parameter identification

`postprocessing/identification/` contains

* `sn_identify.py`: homogeneous S–N life F_p(r) and identification of *p*, η_T and κ from the curves at
  R = −1 and 0. With `sn_experimental_data.csv` it reproduces the parameters of the article (Table 2).
* `paris_identify.py`: identification from reference crack-growth calculations of a CT specimen
  (η_T = 1, R = −1, trial values of *p*) and the measured curves at two stress ratios. The reference
  calculations are run by the user; `paris_experimental_data.csv` contains the measured rates of the
  article. The reference crack histories used in the article are available from the corresponding author
  on request.

```bat
cd postprocessing\identification
python sn_identify.py
python paris_identify.py windows
python paris_identify.py identify <material> <folder with reference crack histories>
```

The reference crack histories are written by `extract_CT_crack_history.py`. The experimental points are
taken from the publications cited in the article.

## 7. Plug-in

`plugin/PF_Fatigue` is an Abaqus/CAE 2023 plug-in that sets up the material, the step, the amplitude, the
output, an optional initial crack and the job for either fatigue law. Installation and use are described
in `plugin/PF-Fatigue_v1.0_User_Guide.pdf`.

## 8. Verification

| Item | Result |
|---|---|
| GC and LE files | differ only in the source term and the flux (line-by-line comparison) |
| `SENT_I_IPC128.inp`, GC and LE (Abaqus 2023, 12 CPUs) | a(N) identical at every output cycle to the runs of the article (GC N = 0–253, LE N = 0–237). These runs used the cycle-by-cycle files before the revision of 12 October 2026, which sets a zero tangent in the clipped branch φ < 0 and keeps the lower bound of *f* in UMATHT (GC file) |
| cycle-domain file vs. the homogeneous closed form (56 single-element S–N cases) | lives 0.18–0.36 % above the closed form; halving ΔN changes them by at most 0.15 % |
| cycle scaling of the cycle-domain law (CT specimens) | lives and rates at different *R* follow from one reference calculation by rescaling the cycle count |
| plug-in kernel (2D cycle by cycle, 2D and 3D cycle domain, error checks) | Abaqus 2023 tests pass |
| plug-in dialog | Check Model, Build and a seam-crack example run in Abaqus/CAE 2023 (User Guide, Section 7) |

## 9. Known issue

`LINK : fatal error LNK1105: cannot close file "standardU.dll"` can appear on Windows when several jobs
compile at the same time or a virus scanner holds the new DLL. The solver has not started in this case;
submit the job again.

## 10. Licence and citation

The package is released under the MIT License (`LICENSE`), with two exceptions. The emblem of Xi'an Jiaotong
University in `plugin/PF_Fatigue/icons/xjtu_logo_*.png` is a trademark of the university; remove or replace these
files when the plug-in is redistributed in modified form. The experimental data in
`postprocessing/identification/*.csv` are taken from Zhang et al. (2016) and Rubio Ruiz et al. (2023) for the
stress-life data and from Zhou et al. (2023) and Ke et al. (2023) for the crack-growth data; please cite these
sources when the data are used. Please cite the article given in `CITATION.cff`.
