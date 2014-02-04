C
C Routine to compute the formal solution of the transfer equation using
C the differencing scheme of Mihalas, Kunasz, and Hummer 1975.
C
C
C                **********************************
C
C If COMPUTE_LAM is true, the approximate LINE diagonal operator is returned.
C This is computed by the method of Schonberg and Hempe (1986 A&A, 163, p151).
C
C Three auxilary arrays are returned. These are
C
C              LAMLINE
C              dZdCHIL
C              dZdETAL
C
C                **********************************
C
C If LINE_BL is true, the routine returns the line EW, the continnum
C intensity.
C
C Three auxilary rays are returned. These are
C
C              JBLANK
C              HBLANK
C          and JINT
C
C JINT must be saved for subsequent calls to FORMSOL. It contains
C the integrated mean intensity of the LINE outside the resonance
C zone. It is required because of our assumption of ``coherent'' electron
C scattering. The mean continuum intensity must also be passed to the
C routine.
C
C                **********************************
C
C 
C
C NB - ETA when passed to this routine must include the electron scattering
C      emissivity.
C    - HQW contains the angular quadrature weights for the H computation
C      at the MID points of the radial mesh.
C    - JQW contains the angular quadrature weights for the J computation
C      on the radial mesh.
C    - The intrinsic line profile is assumed to be depth independent
C      (easily modified but remeber normalization of LFQW).
C
	SUBROUTINE FORMSOL(ETA,CHI,ESEC,CHIL,ETAL,V,SIGMA,R,P,
	1                  JBAR,ZNET,
	1                  LAMLINE,dZdCHIL,dZdETAL,COMPUTE_LAM,
	1                  JCONT,JBLANK,HBLANK,JINT,
	1                  EW,CONT_INT,LINE_BL,FULL_ES,
	1                  JQW,HQW,
	1                  PF,PROF,LFQW,WERFC,FL,DIF,DBB,IC,AMASS,
	1                  THK_LINE,THK_CONT,NLF,NC,NP,ND,METHOD)
C
	IMPLICIT NONE
C
C Altered 28-Oct-1996 - Bug fix: COS corrected back to ACOS in TOR expression.
C Altered 24-May-1996 - SCRTEMP common block removed.
C                       Dynamic allocation now used for temporary arrays.
C                       Generical calls for EXP, COS
C Altered 16-May-1989 - FULL_ES option installed. This option implies
C                       that photons scattered by an electron are not
C                       absorbed by the line.
C Altered 12-May-1989 - Thick boundary condition for EW computaion corrected.
C Altered 03-May-1989 - Base routine was [jdh.lambda]NEWFORMSOL.FOR
C                       Variation of approximate Lambda operator with KI
C                       included - from [jdh.lambda]APPROXNEWT. Bug-fix
C                       Wasn't dividing IBOUND defintion by CHI(1).
C                       Operator calculation installed as an option.
C Altered       -1989 - Blanketing and EW section included.
C Created 31-Jan-1989 - Based on CMFJBAR (and LAMDIAG)
C
C
	INTEGER NLF,NC,NP,ND
	REAL*8 ETA(ND),CHI(ND),ESEC(ND),CHIL(ND),ETAL(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND),P(NP)
	REAL*8 JBAR(ND),ZNET(ND)
	REAL*8 JCONT(ND),JBLANK(ND),HBLANK(ND),JINT(ND)
	REAL*8 LAMLINE(ND),dZdCHIL(ND),dZdETAL(ND)
	REAL*8 JQW(ND,NP),HQW(ND,NP)
	REAL*8 PF(NLF),PROF(NLF),LFQW(NLF),WERFC(NLF)
	REAL*8 DBB,IC,AMASS,FL,EW,CONT_INT
	LOGICAL DIF,THK_CONT,THK_LINE,LINE_BL,FULL_ES,COMPUTE_LAM
	CHARACTER*(*) METHOD
C
	REAL*8 TA(ND),TB(ND),TC(ND),AV(ND),CV(ND),DTAU(ND),Z(ND)
	REAL*8 TCHI(ND),XM(ND),SOURCE(ND),U(ND),VB(ND),VC(ND)
	REAL*8 GB(ND),H(ND),GAM(ND),GAMH(ND),Q(ND),QH(ND)
	REAL*8 JEX_SCAT(ND),BAR_WTS(ND),DEPTH_PRO(ND)
	REAL*8 HCONT(ND),JINT_PREV(ND),HINT(ND)
	REAL*8 AVCONT(ND),CVCONT(ND),dCHIdR(ND)
	REAL*8 AVM1(ND),CVM1(ND),dJBARdKI(ND),DIAG(ND)
	REAL*8 VSRCE(ND),dUdS(ND),dVdS_UP(ND),dVdS_LOW(ND),dJINTdS(ND)
	REAL*8 VKI(ND),RKB(ND),dUdKI(ND),dVdKI_UP(ND),dVdKI_LOW(ND)
C
C Local variables.
C
	INTEGER I,LS,ML,NI,NIEXT
	REAL*8 OLDCHI,T1,T2,DBC,DNU,TOR,IBOUND,WERF_EXP,SCALE
C
C DNU is the frequency bandwidth over which we integrate to obtain
C JINT and is determined by the maximum expansion velocity of the
C atmosphere. Because of the defintion of GAMMA, we define DNU with out
C the factor FL.
C
	DNU=3.33564D-06*V(1)*2.0D0
C
C JINT contains the hand integration of the line from the redmost
C side of the lineprofile (vred) to vred+DNU. If JINT(ND) is zero,
C we assume that this is the first iteration, and set it equal to
C the integral of the continuum intensity.
C
	IF(JINT(ND) .EQ. 0.0D0)THEN
	  DO I=1,ND
	    JINT_PREV(I)=JCONT(I)*DNU
	  END DO
	ELSE
	  DO I=1,ND
	    JINT_PREV(I)=JINT(I)
	  END DO
	END IF
C
C Zero vectors which are to accumulate the integrations.
C
	JBAR(1:ND)=0.0D0
	IF(COMPUTE_LAM)THEN
	  dJBARdKI(1:ND)=0.0D0
	  LAMLINE(1:ND)=0.0D0
	END IF
C
C Vectors are used to compute the Line EW and the blanketed
C MEAN and FLUX intensities.
C
	IF(LINE_BL)THEN
	  JBLANK(1:ND)=0.0D0
	  HBLANK(1:ND)=0.0D0
	  HCONT(1:ND)=0.0D0
	  JINT(1:ND)=0.0D0
	  dJINTdS(1:ND)=0.0D0
	  HINT(1:ND)=0.0D0
	END IF
C
C 
C
C Enter loop to perform integration along each ray.
C
	DO LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)NI=ND
C
C NIEXT is used to compute one extra value of TCHI so that a more accurate
C dCHIdR can be computed.
C
	  NIEXT=NI+1
	  IF(NI .EQ. ND)NIEXT=NI
C
C Zero AV, CV, dUd... and dVd... vectors.
C
	  AV(1:NI)=0.0D0
	  CV(1:NI)=0.0D0
C
	  IF(COMPUTE_LAM)THEN
	    dUdS(1:NI)=0.0D0
	    dVdS_UP(1:NI)=0.0D0
	    dVdS_LOW(1:NI)=0.0D0
	    dUdKI(1:NI)=0.0D0
	    dVdKI_UP(1:NI)=0.0D0
	    dVdKI_LOW(1:NI)=0.0D0
	  END IF
C
	  CALL ZALONGP(R,Z,P(LS),NI)
	  CALL GAMMA(GAM,GAMH,SIGMA,Z,R,V,ND,NI)
C
C Determine boundary condition for continuum intensity.
C
	  IF(THK_CONT)THEN
	    IF(P(LS) .GT. 0.0D0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796D0-ACOS(P(LS)/R(1)))/P(LS)
	    ELSE
	      TOR=CHI(1)*R(1)
	    END IF
	    IBOUND=ETA(1)*(1.0D0-EXP(-TOR))/CHI(1)
	  ELSE
	    TOR=0.0D0
	    IBOUND=0.0D0
	  END IF
C
C 
C
C Perform integration for each frequency in turn.
C
	  OLDCHI=CHI(NI)
	  DO ML=1,NLF
	    T1=PROF(ML)
	    DO I=1,NIEXT
	      DEPTH_PRO(I)=T1	  	  !Used for line profile.
	      TCHI(I)=CHI(I)+CHIL(I)*T1
	      SOURCE(I)=(ETA(I)+ETAL(I)*T1)/TCHI(I)
	    END DO
	    CALL QKIM(Q,QH,GAM,GAMH,TCHI,PF,ML,NI,NLF)
	    IF(DIF .AND. LS .LE. NC)THEN
	      DBC=DBB*DSQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/TCHI(ND)
	1            *(1.0D0+Q(NI)*(1.0D0-TCHI(NI)/OLDCHI))
	    END IF
	    IF(METHOD .EQ. 'ZERO')THEN
	      CALL TAU(DTAU,TCHI,Z,NI)
	    ELSE
	      CALL DERIVCHI(dCHIdR,TCHI,R,NIEXT,METHOD)
	      CALL NORDTAU(DTAU,TCHI,Z,R,dCHIdR,NI)
	    END IF
	    CALL TUVGHD(TA,TB,TC,U,VB,VC,GB,H,Q,QH,DTAU,DIF,LS,NC,NI)
	    CALL XVECD(DTAU,SOURCE,XM,DIF,DBC,IC,LS,NC,ND,NI)
C
C************************************************************************
C   To compute EW need consistency between CONTINUUM and the PRESENT
C   calculation, otherwise we will get spurious equivalent widths.
C   The cases the precice details of the boundary condition are
C   important will be unphysical anyway. For any boundary condition,
C   we expect the flux to be dominated by emission from depth.
C   Important to have boundary condition since strong lines (e.g
C   HeII Lyman Alpha will give a divergent Source function beacuse of the
C   artificial boundary condition TAU=0
C
C   Also need consistency between HCONT and HBLANK.
C
C************************************************************************
C
C If THK_LINE is TRUE, we adopt a "SOBOLEV like" approximation for
C the incident intensity at the outer boundary. The line (pure) incident
C intensity is both angle and frequency independent. The reason for the
C ratio CHI/TCHI is related to the the first order equation for U.
C
C Note that WERFC=-0.5*ERFC where ERFC is the complementary error
C function. We use WERF_EXP because the value for ML=NLF is
C required in the line blanketing section.
C
C	I(incident)=ETAL(1)/CHIL(1)*(1.0D0-WERF_EXP)+IBOUND*WERF_EXP
C
	    IF(THK_LINE)THEN
	      WERF_EXP=EXP(1.0D-15*CHIL(1)/FL/GAM(1)*WERFC(ML))
	      XM(1)=(CHI(1)/TCHI(1)*WERF_EXP-1.0D0)*ETAL(1)/CHIL(1)
	1              -IBOUND*CHI(1)/TCHI(1)*WERF_EXP
	    ELSE
	      XM(1)=-IBOUND
	    END IF
C
C Update AV matrix.
C
	    AV(1)=XM(1)+U(1)*AV(1)
	    DO I=2,NI-1
	      AV(I)=XM(I)+U(I)*AV(I)-VB(I)*CV(I-1)-VC(I)*CV(I)
	    END DO
	    AV(NI)=XM(NI)+U(NI)*AV(NI)-VB(NI)*CV(NI-1)
C
C Save diagonal of TRIDIAGONAL matrix for computation of
C approximate NEWTON-RAPSHON vectors.
C
	    DO I=1,NI
	      DIAG(I)=TB(I)
	    END DO
C
C Solve for the radiation field along ray for this frequency.
C
	    CALL THOMAS(TA,TB,TC,AV,NI,1)
C
C Update C matrix.
C
	    DO I=1,NI-1
	      CV(I)=GB(I)*(AV(I)-AV(I+1))+H(I)*CV(I)
	    END DO
C
C Compute quadrature weights for this frequancy to convert from
C u to JBAR.
C
	    DO I=1,NI
	      BAR_WTS(I)=JQW(I,LS)*LFQW(ML)*DEPTH_PRO(I)
	    END DO
C
C Update JBAR.
C
	    CALL MULTVEC(JBAR,JBAR,BAR_WTS,AV,NI)
C
C 
	    IF(COMPUTE_LAM)THEN
C
C This section of the routine computes the approximate diagonal
C operator. This calculation must be performed before THOMAS is
C called.
C
C Compute the Opacity and the Source variation vectors (diagonal
C operators).
C
	      T1=1.0D-15/FL/GAM(1)*WERFC(ML)
	      CALL VARLAMKI(VKI,VSRCE,RKB,
	1            DEPTH_PRO,DTAU,Z,Q,QH,AV,AVM1,CVM1,
	1            TCHI,SOURCE,ETAL,CHIL,CHI,T1,
	1            THK_LINE,LS,NC,NI)
C
C Update dUdS matrix.
C
	      dUdS(1)=VSRCE(1)+U(1)*dUdS(1)
	      DO I=2,NI-1
	        dUdS(I)=VSRCE(I)+U(I)*dUdS(I)-
	1              VB(I)*dVdS_UP(I-1)-VC(I)*dVdS_LOW(I)
	      END DO
	      dUdS(NI)=VSRCE(NI)+U(NI)*dUdS(NI)-VB(NI)*dVds_UP(NI-1)
C
C Solve for the radiation field along ray for this frequency.
C
	      DO I=1,NI
	        dUdS(I)=dUdS(I)/DIAG(I)
	      END DO
C
C Update dV/dS vectors. Note that V(I) dpends on both S(I) and S(I+1).
C
	      DO I=1,NI-1
	        dVds_UP(I)=-GB(I)*dUdS(I+1)+H(I)*dVds_UP(I)
	        dVds_LOW(I)=GB(I)*dUdS(I)+H(I)*dVds_LOW(I)
	      END DO
C
C 
C
C Update dUdKI matrix.
C
	      dUdKI(1)=VKI(1)+U(1)*dUdKI(1)
	      DO I=2,NI-1
	        dUdKI(I)=VKI(I)+U(I)*dUdKI(I)-
	1               VB(I)*dVdKI_UP(I-1)-VC(I)*dVdKI_LOW(I)
	      END DO
	      dUdKI(NI)=VKI(NI)+U(NI)*dUdKI(NI)-VB(NI)*dVdKI_UP(NI-1)
C
C Solve for the radiation field along ray for this frequency.
C
	      DO I=1,NI
	        dUdKI(I)=dUdKI(I)/DIAG(I)
	      END DO
C
C Update dV/dKI vectors. Note that V(I) depends on both CHIL(I)
c and CHIL(I+1).
C
	      DO I=1,NI-1
	        dVdKI_UP(I)=-GB(I)*dUdKI(I+1)+H(I)*dVdKI_UP(I)+RKB(I)
	        dVdKI_LOW(I)=GB(I)*dUdKI(I)+H(I)*dVdKI_LOW(I)+RKB(I)
	      END DO
C
	      CALL MULTVEC(LAMLINE,LAMLINE,BAR_WTS,dUdS,NI)
	      CALL MULTVEC(dJBARdKI,dJBARdKI,BAR_WTS,dUdKI,NI)
 
	    END IF
C
C 
C
C Save AV (mean intensity lke variable) and CV (mean flux like variable)
C for next frequency integration.
C	
	    DO I=1,NI
	      AVM1(I)=AV(I)
	      CVM1(I)=CV(I)
	    END DO
C
C Update HBLANK and JBLANK.
C
	    IF(LINE_BL)THEN
	      DO I=1,NI
	        JBLANK(I)=JBLANK(I)+AV(I)*JQW(I,LS)*LFQW(ML)
	      END DO
	      DO I=1,NI-1
	        HBLANK(I)=HBLANK(I)+CV(I)*HQW(I,LS)*LFQW(ML)
	      END DO
	      IF(ML .EQ. 1)THEN
	        DO I=1,NI-1
	          HCONT(I)=HCONT(I)+CV(I)*HQW(I,LS)
	          CVCONT(I)=CV(I)
	          AVCONT(I)=AV(I)
	        END DO
	        AVCONT(NI)=AV(NI)
	      END IF
	    END IF
C
	    OLDCHI=TCHI(NI)
	  END DO  		!End DO ML
C
C 
C
	  IF(LINE_BL)THEN
C
	    IF(LS .LE. NC)THEN
	      DBC=DBB*DNU*DSQRT(R(ND)*R(ND)-P(LS)*P(LS))
	1       /R(ND)/CHI(ND)
	    END IF
C
C Compute the total opacity. Store opacity from previous frequency.
C Compute the total source function. We will assume that for the
C the propogation of the line photons that the electron scattering
C is coherent. Since ETA contains a continuum scattering term, we
C first need to subtract this out since it will automatically
C be included in the computations.
C
	    DO I=1,NI
	      SOURCE(I)=( DNU*ETA(I)+ESEC(I)*(JINT_PREV(I)
	1        -JCONT(I)*DNU) )/CHI(I)
	    END DO
C
C Compute Z for this imapct parameter
C
	    CALL ZALONGP(R,Z,P(LS),NI)
	    CALL DERIVCHI(dCHIdR,CHI,R,NIEXT,METHOD)
	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
C
	    DO I=1,NI-1
	      Q(I)=GAM(I)*( AV(I)-AVCONT(I) )/CHI(I)
	      QH(I)=2.0D0*GAMH(I)*( CV(I)-CVCONT(I) )/(CHI(I)+CHI(I+1))
	    END DO
	    Q(NI)=GAM(NI)*( AV(NI)-AVCONT(NI) )/CHI(NI)
C
C For the line, it is very difficult to estimate the incident
C line intensity once where outside the core. We thus assume that
C it is given by the continuum intensity == IBOUND.
C
	    XM(1)=-Q(1)-IBOUND*DNU
	    IF(THK_LINE)XM(1)=XM(1)+GAM(1)*
	1                 (ETAL(1)/CHIL(1)-IBOUND)*(1.0D0-WERF_EXP)
	    VSRCE(1)=0.0D0
	    TA(1)=0.0D0
	    TC(1)=1.0/DTAU(1)
	    TB(1)=-1.0D0-TC(1)
	    DO I=2,NI-1
	      TA(I)=TC(I-1)
	      TC(I)=1.0D0/DTAU(I)
	      TB(I)=-0.5D0*(DTAU(I-1)+DTAU(I))-TA(I)-TC(I)
	      XM(I)=-0.5D0*( SOURCE(I)+Q(I) )*( DTAU(I-1)+DTAU(I) )
	1               -QH(I)+QH(I-1)
	      VSRCE(I)=-0.5D0*( DTAU(I-1)+DTAU(I) )
	    END DO
C
	    IF(LS .LE. NC .AND. DIF)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=TC(NI-1)
	      XM(NI)=DBC			!DNU include in definition of DBC.
	      VSRCE(NI)=0.0D0
	    ELSE IF(LS .GT. NC)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-TA(NI)+DTAU(NI-1)*0.5D0
	      XM(NI)=0.5D0*DTAU(NI-1)*(SOURCE(NI)+Q(NI))-QH(NI-1)
	      VSRCE(NI)=0.5D0*DTAU(NI-1)
	    ELSE
	      TA(NI)=-TC(NI-1)
	      TB(NI)=1-TA(NI)
	      XM(NI)=IC*DNU+Q(NI)
	      VSRCE(NI)=0.0D0
	    END IF
	    TC(NI)=0.0D0
C
C Update dJINT/dS. Since dU/dS=1/TB we must do this before we call THOMAS.
C
	    DO I=1,NI
	      dJINTdS(I)=dJINTdS(I)+VSRCE(I)*JQW(I,LS)/TB(I)
	    END DO
C
C Solve the tridiagonal system of equations.
C
	    CALL THOMAS(TA,TB,TC,XM,NI,1)
C
C Compute the flux
C
	    DO I=1,NI-1
	      CV(I)=(XM(I+1)-XM(I))/DTAU(I)+QH(I)
	    END DO
C
C Update the J and H moments.
C
	    DO I=1,NI
	      JINT(I)=JINT(I)+JQW(I,LS)*XM(I)
	    END DO
	    DO I=1,NI-1
	      HINT(I)=HINT(I)+HQW(I,LS)*CV(I)
	    END DO
C
	  END IF		!End if(LINE_BL)
C 
	END DO			!End do LS
C
	IF(LINE_BL)THEN
C
C Scale JBLANK and HBLANK to allow for the fact that the
C LFQW dont necessarily add to v*(PF(1)-PF(NLF)) (since LFQW is
C normalized so that integral over the line profile is unity).
C
	  SCALE=0.0D0
	  DO ML=1,NLF
	    SCALE=SCALE+LFQW(ML)
	  END DO
	  SCALE=1.0D+15*FL*(PF(1)-PF(NLF))/SCALE
C
	  DO I=1,ND
	    JBLANK(I)=JBLANK(I)*SCALE
	    HBLANK(I)=HBLANK(I)*SCALE
	  END DO
	  IF(FULL_ES)THEN
	    T1=1.0D-15/FL
	    T2=PF(1)-PF(NLF)
	    DO I=1,ND
	      JEX_SCAT(I)=JBLANK(I)*T1-JCONT(I)*T2
	    END DO
	  ELSE
	    JEX_SCAT(1:ND)=0.0D0
	  END IF
C
C Evaluate JBLANK which is deined by Int{Jv dv} over the WHOLE line.
C HBLANK is similarly defined. Note that we need to keep JINT since this
C is used in the electron scattering integral. The frequency factor of
C 10^15 is included in the JBLANK (and HBLANK) definition.
C
	  T1=1.0D+15*FL
	  DO I=1,ND
	    JBLANK(I)=JBLANK(I)+JINT(I)*T1
	    HBLANK(I)=HBLANK(I)+HINT(I)*T1
	  END DO
C
C Update the estimate of JINT using an approximate LAMBDA operator.
C JINT currently contains the FORMAL solution, whilst JINT_PREV
C contains JINT from the previous iteration. NB - JINT is only used
C to compute the electron scattering source function the next time
C through. Thus we add the "Resonance Zone" contibution to it. JEX_SCAT
C is non-zero only if FULL_ES is true.
C
	  DO I=1,ND
	    JINT(I)=JINT(I)+JEX_SCAT(I)
	    T1=dJINTdS(I)*ESEC(I)/CHI(I)
	    JINT(I)=(JINT(I)-T1*JINT_PREV(I))/(1.0D0-T1)
	  END DO
C
C Evaluate the line EW. The units are Angstroms. Also evaluate
C the continuum intensity ( Jys/kpc/kpc ). Note that H is
C defined midway between R(1) and R(2).
C
	  T1=( (PF(1)-PF(NLF))+DNU )*FL*1.0D+15
	  EW=2.99794D-12*( HBLANK(1)-HCONT(1)*T1 )/HCONT(1)/FL/FL
	  CONT_INT=13.19868D0*HCONT(1)*( (R(1)+R(2))**2 )/4.0D0
C
C Change HBLANK and JBLANK to be the int{across line -Jc} (i.e
C integral of J or H above the continuum). Thus, for a weak line,
C HBLANK and JBLANK should be zero.
C
	  DO I=1,ND
	    HBLANK(I)=HBLANK(I)-HCONT(I)*T1
	    JBLANK(I)=JBLANK(I)-JCONT(I)*T1
	  END DO
	END IF
C
C 
C
C Evaluate Z, dZdCHIL and dZdETAL from JBAR, dJBARdS (=LAMLINE) and
C dJBARdKI.
C
	DO I=1,ND
	  ZNET(I)=1.0D0-JBAR(I)*CHIL(I)/ETAL(I)
	  dZdETAL(I)=( JBAR(I)*CHIL(I)/ETAL(I)-LAMLINE(I) )/ETAL(I)
	  dZdCHIL(I)=LAMLINE(I)/CHIL(I) -
	1               ( JBAR(I)+dJBARdKI(I)*CHIL(I) )/ETAL(I)
	END DO
C
	RETURN
	END
