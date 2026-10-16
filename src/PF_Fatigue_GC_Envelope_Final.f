C=======================================================================
C PF_Fatigue_GC_Envelope_Final.f
C Cycle-domain (envelope) phase-field fatigue implementation for Abaqus
C using UMAT + UMATHT.
C
C Temperature DOF is the phase field phi.
C Mechanical response is explicitly staggered: UMAT uses TEMP (phi_n).
C The phase-field equation is the gradient-consistent (GC) form:
C
C   div[f grad(phi)] + 2(1-phi)H/(Gc*l) - f*phi/l^2 = 0
C
C UMATHT defines q = -f grad(phi), so grad(f).grad(phi) is retained by
C the finite-element weak form without explicitly calculating grad(f).
C
C Cycle-domain fatigue law:
C
C   alpha_T0 = Gc/(12*l)
C
C   Amax = g(phi_n)*psi+ / alpha_T0
C
C   Aeq = Amax * ((1-R)/2)^(2*kappa)
C
C   Abar_(N+dN) = Abar_N + dN*Aeq^p
C
C   alpha_T = eta_T*alpha_T0
C
C   f(Abar) = 1,                                  Abar <= eta_T
C           = (2*eta_T/(Abar+eta_T))^2,          Abar >  eta_T
C
C Abar, Amax, Aeq and eta_T are dimensionless.
C DTIME is interpreted as the cycle increment dN in the fatigue step.
C The external Abaqus model is responsible for supplying the peak-load
C envelope state. This subroutine does not detect loading/unloading.
C IMPORTANT: for the envelope law, do not resolve a full triangular cycle
C and simultaneously interpret DTIME as dN. A resolved triangular waveform
C integrates the instantaneous Aeq**p over the cycle and is a different law.
C
C No endurance-limit / fatigue-limit parameter is included.
C
C Mechanical user constants:
C   PROPS(1) = E
C   PROPS(2) = nu
C   PROPS(3) = Gc
C   PROPS(4) = l
C   PROPS(5) = eta_T
C   PROPS(6) = energy split: 0=None, 1=Amor, 2=Miehe, 3=No-tension
C   PROPS(7) = p
C   PROPS(8) = R
C   PROPS(9) = kappa
C
C Residual stiffness is fixed internally:
C   k_res = 1.0D-6
C
C State variables:
C   SDV1 = Abar, dimensionless accumulated fatigue history
C   SDV2 = H, mechanical history field
C   SDV3 = Aeq, current dimensionless equivalent fatigue driver
C   SDV4 = f(Abar), Carrara-type asymptotic degradation
C
C Elasticity, energy split, history field, stress update, phase-field
C source and UMATHT are the same as in PF_Fatigue_Exact_EnergySplit.f;
C only the fatigue law, the user constants and SDV3 differ.
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
      DOUBLE PRECISION E,NU,GC,ELL,ETA,KRES,PEXP,RR,KAPPA
      DOUBLE PRECISION AT0,PHIO,PHIN,GPHI,PSI0,HNEW
      DOUBLE PRECISION PSIP,PSIN
      DOUBLE PRECISION ABAROLD,ABARNEW,AMAX,AEQ,RFAC,FNEW
      INTEGER I,J,ISPLIT
C
      IF (NPROPS .LT. 9) THEN
        WRITE(*,*) '*** PF ENVELOPE: NPROPS must be at least 9.'
        CALL XIT
      END IF
      IF (NSTATV .LT. 4) THEN
        WRITE(*,*) '*** PF ENVELOPE: NSTATV must be at least 4.'
        CALL XIT
      END IF
C
      E      =PROPS(1)
      NU     =PROPS(2)
      GC     =PROPS(3)
      ELL    =PROPS(4)
      ETA    =PROPS(5)
      ISPLIT =NINT(PROPS(6))
      PEXP   =PROPS(7)
      RR     =PROPS(8)
      KAPPA  =PROPS(9)
      KRES   =1.0D-6
C
      IF (GC .LE. 0.0D0 .OR. ELL .LE. 0.0D0 .OR.
     1    ETA .LE. 0.0D0 .OR. PEXP .LE. 0.0D0) THEN
        WRITE(*,*) '*** PF ENVELOPE: Gc, l, eta_T and p must be > 0.'
        CALL XIT
      END IF
      IF (RR .GE. 1.0D0) THEN
        WRITE(*,*) '*** PF ENVELOPE: stress ratio R must be < 1.'
        CALL XIT
      END IF
C
C----- AT2 reference fatigue energy scale ------------------------------
C     alpha_T0 = Gc/(12*l)
C     alpha_T  = eta_T*alpha_T0 is the corresponding dimensional
C     fatigue threshold, while Abar and eta_T are used in dimensionless
C     form in the fatigue accumulation/degradation equations.
      AT0=GC/(12.0D0*ELL)
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
C----- Cycle-domain fatigue accumulation -------------------------------
C     alpha_max = g(phi_n)*psi+
C     Amax      = alpha_max/alpha_T0
C     Aeq       = Amax*((1-R)/2)^(2*kappa)
C     Abar     += dN*Aeq^p, with dN = DTIME
C
C     The external model supplies the peak-load envelope state.
C     No positive-variation detection of alpha is used here.
      ABAROLD=STATEV(1)
      AMAX=GPHI*PSIP/AT0
      RFAC=((1.0D0-RR)/2.0D0)**(2.0D0*KAPPA)
      AEQ=AMAX*RFAC
      IF (AEQ .LT. 0.0D0) AEQ=0.0D0
C
      IF (DTIME .GT. 0.0D0) THEN
        ABARNEW=ABAROLD+DTIME*AEQ**PEXP
      ELSE
        ABARNEW=ABAROLD
      END IF
C
C----- Dimensionless Carrara asymptotic degradation -------------------
C     eta_T = alpha_T/alpha_T0.
      IF (ABARNEW .LE. ETA) THEN
        FNEW=1.0D0
      ELSE
        FNEW=(2.0D0*ETA/(ABARNEW+ETA))**2
      END IF
      IF (FNEW .LT. 1.0D-12) FNEW=1.0D-12
C
C----- GC phase-field local source; diffusion is in UMATHT ------------
      PHIN=TEMP+DTEMP
      IF (PHIN .LT. 0.0D0) THEN
C       Clipped branch: phi=max(TEMP+DTEMP,0) is locally constant.
        PHIN=0.0D0
        RPL=2.0D0*(1.0D0-PHIN)*HNEW/(GC*ELL)
     1      -FNEW*PHIN/(ELL*ELL)
        DRPLDT=0.0D0
      ELSE
C       At PHIN=0 use the right derivative of the admissible phi>=0 branch.
        RPL=2.0D0*(1.0D0-PHIN)*HNEW/(GC*ELL)
     1      -FNEW*PHIN/(ELL*ELL)
        DRPLDT=-2.0D0*HNEW/(GC*ELL)-FNEW/(ELL*ELL)
      END IF
C
C     Mechanical field is intentionally lagged with respect to phi.
      DO I=1,NTENS
        DDSDDT(I)=0.0D0
        DRPLDE(I)=0.0D0
      END DO
C
C----- State variables -------------------------------------------------
      STATEV(1)=ABARNEW
      STATEV(2)=HNEW
      STATEV(3)=AEQ
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
C UMATHT: gradient-consistent phase-field flux q = -f grad(phi)
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
      DOUBLE PRECISION FCOND
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
      FCOND=1.0D0
      IF (NSTATV .GE. 4) THEN
C       STATEV(4)=0 denotes the uninitialized/default state.
C       Once f is active, retain its lower bound instead of reverting to 1
C       when f is exactly 1.0D-12.
        IF (STATEV(4) .GT. 0.0D0) THEN
          FCOND=MAX(STATEV(4),1.0D-12)
        END IF
      END IF
C
      DO I=1,NTGRD
        FLUX(I)=-FCOND*DTEMDX(I)
        DFDG(I,I)=-FCOND
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
