C
C Routine to compute the formal solution of the transfer equation using
C the second order differencing scheme of Hamman (1981)(see also
C Schonburg and Hempe 1986).
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
C If LINE_BL is true, the routine returns the line EW and the observed
C continnum intensity.
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
C Created 31-Jan-1989 - Based on CMFJBAR (and LAMDIAG)
C Altered       -1989 - Blanketing and EW section included.
C Altered 03-May-1989 - Base routine was [jdh.lambda]HAMFORMSOL.FOR
C                       Variation of approximate Lambda operator with KI
C                       included. Bug-fix - wasn't dividing IBOUND defintion
C                       by CHI(1). Operator calculation installed as an option.
C Altered 12-May-1989 - Tick boundary condition for EW section fixed.
C Altered 16-May-1989 - FULL_ES option installed. This option implies
C                       that photons scattered by an electron are not
C                       absorbed by the line.
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
C
	SUBROUTINE HAM_FORMSOL(ETA,CHI,ESEC,CHIL,ETAL,V,SIGMA,R,P,
	1                  JBAR,ZNET,
	1                  LAMLINE,dZdCHIL,dZdETAL,COMPUTE_LAM,
	1                  JCONT,JBLANK,HBLANK,JINT,
	1                  EW,CONT_INT,LINE_BL,FULL_ES,
	1                  JQW,HQW,
	1                  PF,PROF,LFQW,WERFC,FL,DIF,DBB,IC,AMASS,
	1                  THK_LINE,THK_CONT,NLF,NC,NP,ND,METHOD)
C
	IMPLICIT NONE
	INTEGER NLF,NC,NP,ND
	REAL*8 ETA(ND),CHI(ND),ESEC(ND),CHIL(ND),ETAL(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND),P(NP)
	REAL*8 JQW(ND,NP),HQW(ND,NP)
	REAL*8 JBAR(ND),ZNET(ND)
	REAL*8 JCONT(ND),JBLANK(ND),HBLANK(ND),JINT(ND)
	REAL*8 LAMLINE(ND),dZdCHIL(ND),dZdETAL(ND)
	REAL*8 PF(NLF),PROF(NLF),LFQW(NLF),WERFC(NLF)
	REAL*8 DBB,IC,AMASS,FL,EW,CONT_INT
	LOGICAL DIF,THK_CONT,THK_LINE,LINE_BL,FULL_ES,COMPUTE_LAM
	CHARACTER*(*) METHOD
C
	INTEGER NV
	PARAMETER (NV=100)
	COMMON /SCRATCH/  PROGDESC,TA,TB,TC,AV,CV,DTAU,Z,TCHI,XM,SOURCE,
	1                   UA,UB,UC,
	1                   VB,VC,GB,H,GAM,GAMH,Q,QH,
	1                   JEX_SCAT,BAR_WTS,DEPTH_PRO,
	1                   AVM1,CVM1,dJBARdKI,DIAG,
	1                   HCONT,JINT_PREV,HINT,
	1                   AVCONT,CVCONT,dCHIdR,
	1                   VSRCE,dUdS,dVdS_UP,dVdS_LOW,dJINTdS,
	1                   VKI,RKB,dUdKI,dVdKI_UP,dVdKI_LOW
C
	REAL*8 PROGDESC
	REAL*8 TA(NV),TB(NV),TC(NV),AV(NV),CV(NV),DTAU(NV),Z(NV)
	REAL*8 UA(NV),UB(NV),UC(NV)
	REAL*8 TCHI(NV),XM(NV),SOURCE(NV),VB(NV),VC(NV)
	REAL*8 GB(NV),H(NV),GAM(NV),GAMH(NV),Q(NV),QH(NV)
	REAL*8 JEX_SCAT(NV),BAR_WTS(NV),DEPTH_PRO(NV)
	REAL*8 AVM1(NV),CVM1(NV),dJBARdKI(NV),DIAG(NV)
	REAL*8 dJINTdS(NV),HCONT(NV),JINT_PREV(NV),HINT(NV)
	REAL*8 AVCONT(NV),CVCONT(NV),dCHIdR
	REAL*8 VSRCE(NV),dUdS(NV),dVdS_UP(NV),dVdS_LOW(NV)
	REAL*8 VKI(NV),RKB(NV),dUdKI(NV),dVdKI_UP(NV),dVdKI_LOW(NV)
C
C Local variables.
C
	INTEGER I,LS,ML,NI,NIEXT
	REAL*8 OLDCHI,T1,T2,DBC,VDBC,DNU,TOR,IBOUND,MID_PRO,SCALE
	REAL*8 WERF_EXP
C
C PROGDESC is a variable use to confirm that the scratch block is not
C being used by some other routine.
C
	PROGDESC=79765231.0D0		!Must be unique
C
	IF(ND .GT. NV)THEN
	  WRITE(6,*)'Error in HAM_FORMSOL - NV not large enough'
	  STOP
	END IF
C
C DNU is the frequency bandwidth over which we integrate to obtain
C JINT and is determined by the maximum expansion velocity of the
C atmosphere. We define DNU with out the factor FL.
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
	CALL DP_ZERO(JBAR,ND)
	IF(COMPUTE_LAM)THEN
	  CALL DP_ZERO(dJBARdKI,ND)
	  CALL DP_ZERO(LAMLINE,ND)
	END IF
C
C Vectors are used to compute the Line EW and the blanketed
C MEAN and FLUX intensities.
C
	IF(LINE_BL)THEN
	  CALL DP_ZERO(JBLANK,ND)
	  CALL DP_ZERO(HBLANK,ND)
	  CALL DP_ZERO(HCONT,ND)
	  CALL DP_ZERO(JINT,ND)
	  CALL DP_ZERO(dJINTdS,ND)
	  CALL DP_ZERO(HINT,ND)
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
C NIEXT is used to compute one extral value of TCHI so that a more accurate
C dCHIdR can be computed.
C
	  NIEXT=NI+1
	  IF(NI .EQ. ND)NIEXT=NI
C
C Zero important vectors.
C
	  CALL DP_ZERO(AV,NI)
	  CALL DP_ZERO(CV,NI)
	  CALL DP_ZERO(AVM1,NI)
	  CALL DP_ZERO(CVM1,NI)
	  IF(COMPUTE_LAM)THEN
	    CALL DP_ZERO(dUdS,NI)
	    CALL DP_ZERO(dVdS_UP,NI)
	    CALL DP_ZERO(dVdS_LOW,NI)
	    CALL DP_ZERO(dUdKI,NI)
	    CALL DP_ZERO(dVdKI_UP,NI)
	    CALL DP_ZERO(dVdKI_LOW,NI)
	  END IF
C
	  CALL ZALONGP(R,Z,P(LS),NI)
	  CALL GAMMA(GAM,GAMH,SIGMA,Z,R,V,ND,NI)
C
C Determine boundary condition for continuum intensity.
C
	  IF(THK_CONT)THEN
	    IF(P(LS) .GT. 0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796-DACOS(P(LS)/R(1)))/P(LS)
	    ELSE
	      TOR=CHI(1)*R(1)
	    END IF
	    IBOUND=ETA(1)*(1.0D0-DEXP(-TOR))/CHI(1)
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
	    IF(ML .EQ. 1)THEN
	      MID_PRO=0.0D0
	    ELSE
	      MID_PRO=0.50D0*( PROF(ML)+PROF(ML-1) )
	    END IF
	    DO I=1,NIEXT
	      DEPTH_PRO(I)=MID_PRO	  	  !Used for dJBARdS
	      TCHI(I)=CHI(I)+CHIL(I)*MID_PRO
	      SOURCE(I)=(ETA(I)+ETAL(I)*MID_PRO)/TCHI(I)
	    END DO
	    CALL QKIM(Q,QH,GAM,GAMH,TCHI,PF,ML,NI,NLF)
	    IF(DIF .AND. LS .LE. NC)THEN
	      T1=CHI(ND)+CHIL(ND)*PROF(ML)
	      DBC=DBB*DSQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/TCHI(ND)
	1            *( 1.0D0+Q(ND)*TCHI(ND)*(1.0D0/T1-1.0D0/OLDCHI) )
C
C Variation with respect to T1 and OLDCHI only.
C
	      IF(ML .EQ. 1)THEN
	        VDBC=0.0D0
	      ELSE
	        VDBC=DBB*Q(ND)*(PROF(ML-1)/OLDCHI/OLDCHI-PROF(ML)/T1/T1)
	1            *DSQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)
	      END IF
	      OLDCHI=T1
	    END IF
	    IF(METHOD .EQ. 'ZERO')THEN
	      CALL TAU(DTAU,TCHI,Z,NI)
	    ELSE
	      CALL DERIVCHI(dCHIdR,TCHI,R,NIEXT,METHOD)
	      CALL NORDTAU(DTAU,TCHI,Z,R,dCHIdR,NI)
	    END IF
	    CALL HAMTUVGH(TA,TB,TC,UA,UB,UC,VB,VC,GB,H,XM,Q,QH,
	1                   SOURCE,DTAU,DBC,IC,DIF,LS,NC,NI,ML)
C
C If THK_LINE is TRUE, we adopt a "SOBOLEV like" approximation for
C the incident intensity at the outer boundary. The line (pure) incident
C intensity is both angle and frequency independent. The reason for the
C ratio CHI/TCHI is related to the the first order equation for U.
C
C Note that WERFC=-0.5*ERFC where ERFC is the complementary error
C function.
C
C	I(incident)=ETAL(I)/CHIL(I)*(1.0D0-WERF_EXP)+IBOUND*WERF_EXP
C
	    IF(THK_LINE .AND. ML .NE. 1)THEN
	      WERF_EXP=DEXP( 1.0D-15*CHIL(1)/FL/GAM(1)*
	1                      0.5D0*(WERFC(ML)+WERFC(ML-1)) )
	      XM(1)=(CHI(1)/TCHI(1)*WERF_EXP-1.0D0)*ETAL(1)/CHIL(1)
	1              -IBOUND*CHI(1)/TCHI(1)*WERF_EXP
	    ELSE
	      XM(1)=-IBOUND
	    END IF
C
C Update SOURCE vector for coupling to the previous frequency.
C As AV is still required to update CV, we store result in XM.
C
	    AV(1)=XM(1)+UB(1)*AVM1(1)+UC(1)*AVM1(2)+VC(1)*CVM1(1)
	    DO I=2,NI-1
	      AV(I)=XM(I)+UA(I)*AVM1(I-1)+UB(I)*AVM1(I)+UC(I)*AVM1(I+1)
	1           +VB(I)*CVM1(I-1)+VC(I)*CVM1(I)
	    END DO
	    AV(NI)=XM(NI)+UA(NI)*AVM1(NI-1)+UB(NI)*AVM1(NI)
	1           +VB(NI)*CVM1(NI-1)
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
C Update CV vector. NB Sign of GB is reveresd to that in CMFJBAR.
C
	    DO I=1,NI-1
	      CV(I)=GB(I)*(AV(I+1)-AV(I)+AVM1(I+1)-AVM1(I))+H(I)*CVM1(I)
	    END DO
C
C Compute quadrature weights for this frequency to convert from
C u to JBAR.
C
	    DO I=1,NI
	      BAR_WTS(I)=JQW(I,LS)*LFQW(ML)*PROF(ML)
	    END DO
C
C Update JBAR.
C
	    CALL MULTVEC(JBAR,JBAR,BAR_WTS,AV,NI)
C
C 
C
C This section of the routine computes the approximate diagonal
C operator. Note that dU/dSL and dU/dKIL are intrisically zero
C for the first frequency.
C
	    IF(COMPUTE_LAM .AND. ML .NE. 1)THEN
C
C Compute the Opacity and the Source variation vectors (diagonal
C operators).
C
	      T1=1.0D-15/FL/GAM(1)*( WERFC(ML)+WERFC(ML-1) )*0.5D0
	      CALL HAM_VARLAM(VKI,VSRCE,RKB,
	1              DEPTH_PRO,DTAU,R,Z,
	1              GB,Q,QH,AV,AVM1,CVM1,
	1              TCHI,SOURCE,ETAL,CHIL,CHI,
	1              IBOUND,T1,THK_LINE,
	1              DBC,DIF,LS,NC,NI)
C
C Update VKI if diffusion approximation. Note that VDBC refers only to
C the derivative of the TCHI*Q*( - ) in the expression for DBC.
C
	      IF(LS .LE. NC .AND. DIF)THEN
	        VKI(NI)=VKI(NI)+VDBC
	      END IF
C
C Update dUdS matrix.
C
	      VSRCE(1)=VSRCE(1)+UB(1)*dUdS(1)+VC(1)*dVdS_LOW(1)
	      DO I=2,NI-1
	        VSRCE(I)=VSRCE(I)+UB(I)*dUdS(I)
	1             +VB(I)*dVdS_UP(I-1)+VC(I)*dVdS_LOW(I)
	      END DO
	      VSRCE(NI)=VSRCE(NI)+UB(NI)*dUdS(NI)
	1             +VB(NI)*dVdS_UP(NI-1)
C
C Solve for the radiation field along ray for this frequency.
C
	      DO I=1,NI
	        VSRCE(I)=VSRCE(I)/DIAG(I)
	      END DO
C
C Update dV/dS vectors. Note that V(I) dpends on both S(I) and S(I+1).
C Also, we have changed the sign of GB form CMFJBAR.
C
	      DO I=1,NI-1
	        dVds_UP(I)=H(I)*dVds_UP(I)
	1                    +GB(I)*( dUdS(I+1)+VSRCE(I+1) )
	        dVds_LOW(I)=H(I)*dVds_LOW(I)
	1                    -GB(I)*( dUdS(I)+VSRCE(I) )
	      END DO
C
C As don't require dUdS at previous frequency anymore, we can update its
C values to that of the current frequency.
C
	      DO I=1,NI
	        dUdS(I)=VSRCE(I)
	      END DO
C 
C
C Update dUdKI matrix.
C
	      VKI(1)=VKI(1)+UB(1)*dUdKI(1)+VC(1)*dVdKI_LOW(1)
	      DO I=2,NI-1
	        VKI(I)=VKI(I)+UB(I)*dUdKI(I)
	1             +VB(I)*dVdKI_LOW(I-1)+VC(I)*dVdKI_UP(I)
	      END DO
	      VKI(NI)=VKI(NI)+UB(NI)*dUdKI(NI)
	1             +VB(NI)*dVdKI_UP(NI-1)
C
C Solve for the radiation field along ray for this frequency.
C
	      DO I=1,NI
	        VKI(I)=VKI(I)/DIAG(I)
	      END DO
C
C Update dV/dKI vectors. Note that V(I) depends on both CHIL(I)
c and CHIL(I+1).
C
	      DO I=1,NI-1
	        dVdKI_UP(I)=H(I)*dVdKI_UP(I)+RKB(I)
	1           +GB(I)*( dUdKI(I+1)+VKI(I+1) )
	        dVdKI_LOW(I)=H(I)*dVdKI_LOW(I)+RKB(I)
	1           -GB(I)*( dUdKI(I)+VKI(I) )
	      END DO
C
C As don't require dUdS at previous frequency anymore, we can update its
C values to that of the current frequency.
C
	      DO I=1,NI
	        dUdKI(I)=VKI(I)
	      END DO
C
C Update diagonal operators.
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
	    IF(LINE_BL)THEN
C
C Update HBLANK and JBLANK.
C
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
	    IF(METHOD .EQ. 'ZERO')THEN
	      CALL TAU(DTAU,CHI,Z,NI)
	    ELSE
	      CALL DERIVCHI(dCHIdR,CHI,R,NIEXT,METHOD)
	      CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
	    END IF
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
C Note that WERF_EXP must be evaluate at NLF - not NLF-0.5 .
C
	    XM(1)=-Q(1)-IBOUND*DNU
	    IF(THK_LINE)THEN
	      WERF_EXP=DEXP( 1.0D-15*CHIL(1)/FL/GAM(1)*WERFC(NLF) )
	      XM(1)=XM(1)+GAM(1)*
	1                 (ETAL(1)/CHIL(1)-IBOUND)*(1.0D0-WERF_EXP)
	    END IF
	    VSRCE(1)=0.0D0
	    TA(1)=0.0D0
	    TC(1)=1./DTAU(1)
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
	    TC(NI)=0.0
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
	    CALL DP_ZERO(JEX_SCAT,ND)
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
	  CONT_INT=13.19868*HCONT(1)*( (R(1)+R(2))**2 )/4.0D0
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
C Evaluate ZNET, dZdCHIL and dZdETAL from JBAR, dJBARdS (=LAMLINE) and
C dJBARdKI.
C
	DO I=1,ND
	  ZNET(I)=1.0D0-JBAR(I)*CHIL(I)/ETAL(I)
	END DO
C
	IF(COMPUTE_LAM)THEN
	  DO I=1,ND
	    dZdETAL(I)=( JBAR(I)*CHIL(I)/ETAL(I)-LAMLINE(I) )/ETAL(I)
	    dZdCHIL(I)=LAMLINE(I)/CHIL(I) -
	1               ( JBAR(I)+dJBARdKI(I)*CHIL(I) )/ETAL(I)
	  END DO
	END IF
C
	IF(PROGDESC .NE. 79765231.0D0)THEN
	  WRITE(6,*)'Error - Scratch Block Corrupted in FORMSOL'
	  STOP
	END IF
C
	RETURN
	END
