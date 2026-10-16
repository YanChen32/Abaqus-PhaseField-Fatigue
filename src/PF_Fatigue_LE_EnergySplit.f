C=======================================================================
C PF_Fatigue_LE_EnergySplit.f
C Minimal cycle-by-cycle phase-field fatigue implementation for Abaqus
C Local-equivalent (LE) comparison version using UMAT + UMATHT.
C
C Revision 2026-10-12: clipped-branch tangent aligned with
C PF_Fatigue_GC_Envelope_Final.f; header comments corrected.
C
C Temperature DOF is the phase field phi.
C Mechanical response is explicitly staggered: UMAT uses TEMP (phi_n).
C The phase-field equation uses the current H and f after mechanics:
C
C   lap(phi) + 2(1-phi)H/(f*Gc*l) - phi/l^2 = 0
C
C UMATHT defines q = -grad(phi); the term grad(ln f).grad(phi) of the
C gradient-consistent (GC) equation is therefore not represented.
C
C Mechanical user constants:
C   PROPS(1) = E
C   PROPS(2) = nu
C   PROPS(3) = Gc
C   PROPS(4) = l
C   PROPS(5) = alpha_T
C   PROPS(6) = residual stiffness k_res
C   PROPS(7) = energy split: 0=None, 1=Amor, 2=Miehe, 3=No-tension
C
C State variables:
C   SDV1 = alpha_bar
C   SDV2 = H
C   SDV3 = alpha at previous converged increment
C   SDV4 = f(alpha_bar), Carrara asymptotic degradation
C
C Isotropic AT2 with cycle-by-cycle fatigue accumulation.
C=======================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),
     4 DFGRD1(3,3)
C
      DOUBLE PRECISION C0(6,6),EPS(6),SIG0(6)
      DOUBLE PRECISION E,NU,GC,ELL,AT,KRES
      DOUBLE PRECISION PHIO,PHIN,GPHI,PSI0,HNEW
      DOUBLE PRECISION PSIP,PSIN
      DOUBLE PRECISION AOLD,ACUR,APREV,ANEW,FNEW,DALPHA
      INTEGER I,J,ISPLIT
C
      IF (NPROPS .LT. 7) THEN
        WRITE(*,*) '*** PF LE: NPROPS must be at least 7.'
        CALL XIT
      END IF
      IF (NSTATV .LT. 4) THEN
        WRITE(*,*) '*** PF LE: NSTATV must be at least 4.'
        CALL XIT
      END IF
C
      E    = PROPS(1)
      NU   = PROPS(2)
      GC   = PROPS(3)
      ELL  = PROPS(4)
      AT   = PROPS(5)
      KRES = PROPS(6)
      ISPLIT=NINT(PROPS(7))
      IF (KRES .LT. 0.0D0) KRES = 0.0D0
C
      IF (GC .LE. 0.0D0 .OR. ELL .LE. 0.0D0 .OR.
     1    AT .LE. 0.0D0) THEN
        WRITE(*,*) '*** PF LE: Gc, l and alpha_T must be > 0.'
        CALL XIT
      END IF
C
C----- Undamaged elastic matrix ----------------------------------------
      CALL PF_ELASTIC(C0,E,NU,NDI,NSHR,NTENS)
C
      DO I=1,6
        EPS(I)=0.0D0
        SIG0(I)=0.0D0
      END DO
      DO I=1,NTENS
        EPS(I)=STRAN(I)+DSTRAN(I)
      END DO
      DO I=1,NTENS
        DO J=1,NTENS
          SIG0(I)=SIG0(I)+C0(I,J)*EPS(J)
        END DO
      END DO
C
C----- Explicit stagger: mechanics uses phase at increment start -------
      PHIO=TEMP
      IF (PHIO .LT. 0.0D0) PHIO=0.0D0
      IF (PHIO .GT. 1.0D0) PHIO=1.0D0
      GPHI=(1.0D0-PHIO)**2+KRES
C
      DO I=1,NTENS
        STRESS(I)=GPHI*SIG0(I)
        DO J=1,NTENS
          DDSDDE(I,J)=GPHI*C0(I,J)
        END DO
      END DO
C
C----- Undegraded elastic energy and history field --------------------
      PSI0=0.0D0
      DO I=1,NTENS
        PSI0=PSI0+0.5D0*SIG0(I)*EPS(I)
      END DO
      IF (PSI0 .LT. 0.0D0) PSI0=0.0D0
C
C----- Energy split: 0=None, 1=Amor, 2=Miehe, 3=No-tension ----------
      CALL PF_ENERGY_SPLIT(EPS,E,NU,NDI,NSHR,NTENS,ISPLIT,
     1                     PSI0,PSIP,PSIN)
C
      HNEW=MAX(STATEV(2),PSIP)
C
C----- Cycle-by-cycle fatigue accumulation ----------------------------
C     alpha = g(phi_n)*psi+
C     alpha_bar += positive variation of alpha only.
      AOLD =STATEV(1)
      APREV=STATEV(3)
      ACUR =GPHI*PSIP
      DALPHA=ACUR-APREV
      IF (DALPHA .GT. 0.0D0 .AND. DTIME .GT. 0.0D0) THEN
        ANEW=AOLD+DALPHA
      ELSE
        ANEW=AOLD
      END IF
C
C----- Carrara asymptotic degradation ---------------------------------
      IF (ANEW .LE. AT) THEN
        FNEW=1.0D0
      ELSE
        FNEW=(2.0D0*AT/(ANEW+AT))**2
      END IF
      IF (FNEW .LT. 1.0D-12) FNEW=1.0D-12
C
C----- LE phase-field local source; diffusion is in UMATHT ------------
      PHIN=TEMP+DTEMP
      IF (PHIN .LT. 0.0D0) THEN
C       Clipped branch: phi=max(TEMP+DTEMP,0) is locally constant.
        PHIN=0.0D0
        RPL=2.0D0*(1.0D0-PHIN)*HNEW/(FNEW*GC*ELL)
     1      -PHIN/(ELL*ELL)
        DRPLDT=0.0D0
      ELSE
C       At PHIN=0 use the right derivative of the admissible phi>=0 branch.
        RPL=2.0D0*(1.0D0-PHIN)*HNEW/(FNEW*GC*ELL)
     1      -PHIN/(ELL*ELL)
        DRPLDT=-2.0D0*HNEW/(FNEW*GC*ELL)-1.0D0/(ELL*ELL)
      END IF
C
C     Mechanical field is intentionally lagged with respect to phi.
      DO I=1,NTENS
        DDSDDT(I)=0.0D0
        DRPLDE(I)=0.0D0
      END DO
C
C----- State variables -------------------------------------------------
      STATEV(1)=ANEW
      STATEV(2)=HNEW
      STATEV(3)=ACUR
      STATEV(4)=FNEW
C
      SSE=GPHI*PSI0
      SPD=0.0D0
      SCD=0.0D0
C
      RETURN
      END
C
C=======================================================================
C UMATHT: constant-coefficient (LE) phase-field flux
C=======================================================================
      SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1 STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2 CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3 NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1 DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)
C
      INTEGER I,J
C
      U=0.0D0
      DUDT=0.0D0
      DO I=1,NTGRD
        DUDG(I)=0.0D0
        DFDT(I)=0.0D0
        DO J=1,NTGRD
          DFDG(I,J)=0.0D0
        END DO
      END DO
C
      DO I=1,NTGRD
        FLUX(I)=-DTEMDX(I)
        DFDG(I,I)=-1.0D0
      END DO
C
      RETURN
      END
C
C=======================================================================
C Isotropic linear-elastic matrix for plane stress, plane strain/
C axisymmetric, and 3D small-strain elements.
C=======================================================================
      SUBROUTINE PF_ELASTIC(C0,E,NU,NDI,NSHR,NTENS)
      INCLUDE 'ABA_PARAM.INC'
      DOUBLE PRECISION C0(6,6),E,NU,LAM,MU,FAC
      INTEGER NDI,NSHR,NTENS,I,J
C
      DO I=1,6
        DO J=1,6
          C0(I,J)=0.0D0
        END DO
      END DO
C
      IF (NDI .EQ. 2 .AND. NTENS .EQ. 3) THEN
C       Plane stress
        FAC=E/(1.0D0-NU*NU)
        C0(1,1)=FAC
        C0(2,2)=FAC
        C0(1,2)=FAC*NU
        C0(2,1)=FAC*NU
        C0(3,3)=E/(2.0D0*(1.0D0+NU))
      ELSE IF (NDI .EQ. 3) THEN
C       3D, plane strain, or axisymmetric
        MU=E/(2.0D0*(1.0D0+NU))
        LAM=E*NU/((1.0D0+NU)*(1.0D0-2.0D0*NU))
        DO I=1,3
          DO J=1,3
            C0(I,J)=LAM
          END DO
          C0(I,I)=LAM+2.0D0*MU
        END DO
        DO I=4,NTENS
          C0(I,I)=MU
        END DO
      ELSE
        WRITE(*,*) '*** PF: unsupported stress state. NDI=',NDI,
     1             ' NTENS=',NTENS
        CALL XIT
      END IF
C
      RETURN
      END
C
C=======================================================================
C Energy split used only for fracture/fatigue driving energy.
C Stress and tangent are degraded as a whole (hybrid formulation).
C
C ISPLIT=0: no split
C ISPLIT=1: Amor volumetric-deviatoric split
C ISPLIT=2: Miehe spectral split
C ISPLIT=3: No-tension split
C=======================================================================
      SUBROUTINE PF_ENERGY_SPLIT(EPS,E,NU,NDI,NSHR,NTENS,
     1                           ISPLIT,PSI0,PSIP,PSIN)
      INCLUDE 'ABA_PARAM.INC'
C
      DOUBLE PRECISION EPS(6),E,NU,PSI0,PSIP,PSIN
      DOUBLE PRECISION EWORK(6),PS(3),EDEV(3)
      DOUBLE PRECISION BK,EG,TRE,EDEVS,TREP,TREN
      DOUBLE PRECISION TRP1,TRN1,TRP2,TRN2
      DOUBLE PRECISION AKK,BKK,CKK
      INTEGER NDI,NSHR,NTENS,ISPLIT,I,NSHR3
C
C----- No split ---------------------------------------------------------
      IF (ISPLIT .EQ. 0) THEN
        PSIP=PSI0
        PSIN=0.0D0
        RETURN
      END IF
C
      IF (ISPLIT .LT. 0 .OR. ISPLIT .GT. 3) THEN
        WRITE(*,*) '*** PF: invalid energy split selector=',ISPLIT
        CALL XIT
      END IF
C
C----- Build a 3D strain vector for principal-strain evaluation --------
C     For plane stress, reconstruct epsilon_33 from sigma_33=0.
      DO I=1,6
        EWORK(I)=0.0D0
      END DO
C
      IF (NDI .EQ. 2 .AND. NTENS .EQ. 3) THEN
        EWORK(1)=EPS(1)
        EWORK(2)=EPS(2)
        EWORK(3)=-NU*(EPS(1)+EPS(2))/(1.0D0-NU)
        EWORK(4)=EPS(3)
        NSHR3=1
      ELSE
        DO I=1,NTENS
          EWORK(I)=EPS(I)
        END DO
        NSHR3=NSHR
      END IF
C
      CALL SPRINC(EWORK,PS,2,3,NSHR3)
C
      BK=E/(1.0D0-2.0D0*NU)/3.0D0
      EG=E/(1.0D0+NU)/2.0D0
C
C----- Amor volumetric-deviatoric split -------------------------------
      IF (ISPLIT .EQ. 1) THEN
        TRE=PS(1)+PS(2)+PS(3)
        EDEV(1)=PS(1)-TRE/3.0D0
        EDEV(2)=PS(2)-TRE/3.0D0
        EDEV(3)=PS(3)-TRE/3.0D0
        EDEVS=EDEV(1)**2+EDEV(2)**2+EDEV(3)**2
        TREP=0.5D0*(TRE+ABS(TRE))
        TREN=0.5D0*(TRE-ABS(TRE))
        PSIP=0.5D0*BK*TREP**2+EG*EDEVS
        PSIN=0.5D0*BK*TREN**2
        RETURN
      END IF
C
C----- Miehe spectral split -------------------------------------------
      IF (ISPLIT .EQ. 2) THEN
        TRP1=(PS(1)+PS(2)+PS(3)
     1       +ABS(PS(1)+PS(2)+PS(3)))/2.0D0
        TRN1=(PS(1)+PS(2)+PS(3)
     1       -ABS(PS(1)+PS(2)+PS(3)))/2.0D0
        TRP2=0.0D0
        TRN2=0.0D0
        DO I=1,3
          TRP2=TRP2+(PS(I)+ABS(PS(I)))**2/4.0D0
          TRN2=TRN2+(PS(I)-ABS(PS(I)))**2/4.0D0
        END DO
        PSIP=NU*EG/(1.0D0-2.0D0*NU)*TRP1**2+EG*TRP2
        PSIN=NU*EG/(1.0D0-2.0D0*NU)*TRN1**2+EG*TRN2
        RETURN
      END IF
C
C----- No-tension split ------------------------------------------------
C     Sort principal strains: PS(1)>=PS(2)>=PS(3).
      AKK=MAX(PS(1),PS(2),PS(3))
      BKK=MIN(PS(1),PS(2),PS(3))
      CKK=PS(1)+PS(2)+PS(3)-AKK-BKK
      PS(1)=AKK
      PS(2)=CKK
      PS(3)=BKK
C
      IF (PS(3) .GT. 0.0D0) THEN
        PSIP=NU*EG/(1.0D0-2.0D0*NU)
     1       *(PS(1)+PS(2)+PS(3))**2
     2       +EG*(PS(1)**2+PS(2)**2+PS(3)**2)
        PSIN=0.0D0
      ELSE IF ((PS(2)+NU*PS(3)) .GT. 0.0D0) THEN
        PSIP=NU*EG/(1.0D0-2.0D0*NU)
     1       *(PS(1)+2.0D0*NU*PS(3)+PS(2))**2
     2       +EG*((PS(1)+NU*PS(3))**2
     3       +(PS(2)+NU*PS(3))**2)
        PSIN=E/2.0D0*PS(3)**2
      ELSE IF (((1.0D0-NU)*PS(1)
     1          +NU*(PS(2)+PS(3))) .GT. 0.0D0) THEN
        PSIP=E/(2.0D0*((1.0D0-NU**2)*(1.0D0-2.0D0*NU)))
     1       *((1.0D0-NU)*PS(1)+NU*PS(2)+NU*PS(3))**2
        PSIN=E/(2.0D0*(1.0D0-NU**2))
     1       *(PS(2)**2+PS(3)**2+2.0D0*NU*PS(2)*PS(3))
      ELSE
        PSIP=0.0D0
        PSIN=NU*EG/(1.0D0-2.0D0*NU)
     1       *(PS(1)+PS(2)+PS(3))**2
     2       +EG*(PS(1)**2+PS(2)**2+PS(3)**2)
      END IF
C
      RETURN
      END
