C
C Routine to compute the Eddington F and G factors. The second order
C differencing scheme of Hamman (1981)(see also Schonburg and Hempe 1986)
C is used.
C
C                **********************************
C                **********************************
C
C NB - ETA when passed to this routine must include the electron scattering
C      emissivity.
C    - HQW, NQW contain the angular quadrature weights for the H, N computation
C      at the MID points of the radial mesh.
C    - JQW, KQW contain the angular quadrature weights for the J, K computation
C      on the radial mesh.
C    - The intrinsic line profile is assumed to be depth independent
C      (easily modified but remember normalization of LFQW).
C
C Created 10-May-1989 Based on FG_COMP and HAM_FORMSOL.
C Altered 12-May-1989 - Thick boundary condition in EW section fixed.
C Altered 02-Jun-1989 - Thick boundary correction corrected. HBC was
C                       being computed incorrectly.
C
	SUBROUTINE FG_HAM(ETA,CHI,ESEC,JCONT,
	1                  CHIL,ETAL,V,SIGMA,R,P,
	1                  JNU,HNU,KNU,NNU,
	1                  JQW,HQW,KQW,NQW,
	1                  IN_HBC,HBC,NBC,
	1                  PF,PROF,LFQW,WERFC,FL,
	1                  EW,CONT_INT,LINE_BL,
	1                  DIF,DBB,IC,METHOD,
	1                  THK_CONT,THK_LINE,NLF,NC,NP,ND)
	IMPLICIT NONE
C
	INTEGER NLF,NC,NP,ND
	REAL*8 ETA(ND),CHI(ND),ESEC(ND),JCONT(ND),CHIL(ND),ETAL(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND),P(NP)
C
	REAL*8 JQW(ND,NP),HQW(ND,NP),KQW(ND,NP),NQW(ND,NP)
	REAL*8 JNU(ND,NLF+1),HNU(ND,NLF+1),KNU(ND,NLF+1),NNU(ND,NLF+1)
	REAL*8 IN_HBC(NLF+1),HBC(NLF+1),NBC(NLF+1)
C
	REAL*8 PF(NLF),PROF(NLF),LFQW(NLF),WERFC(NLF)
	REAL*8 DBB,IC,FL,EW,CONT_INT
	LOGICAL DIF,LINE_BL,THK_CONT,THK_LINE
	CHARACTER*(*) METHOD
C
	INTEGER NV
	PARAMETER (NV=100)
	COMMON /SCRATCH/  PROGDESC,TA,TB,TC,AV,CV,DTAU,Z,TCHI,XM,SOURCE,
	1                   UA,UB,UC,
	1                   VB,VC,GB,H,GAM,GAMH,Q,QH,
	1                   AVM1,CVM1,JPREV,
	1                   AVCONT,CVCONT,dCHIdR
C
	REAL*8 PROGDESC
	REAL*8 TA(NV),TB(NV),TC(NV),AV(NV),CV(NV),DTAU(NV),Z(NV)
	REAL*8 TCHI(NV),XM(NV),SOURCE(NV)
	REAL*8 UA(NV),UB(NV),UC(NV),VB(NV),VC(NV)
	REAL*8 GB(NV),H(NV),GAM(NV),GAMH(NV),Q(NV),QH(NV)
	REAL*8 AVM1(NV),CVM1(NV),JPREV(NV)
	REAL*8 AVCONT(NV),CVCONT(NV),dCHIdR(NV)
C
C Local variables.
C
	INTEGER I,J,LS,ML,NI,NIEXT
	REAL*8 OLDCHI,T1,DBC,DNU,TOR,MID_PRO,SCALE,HBLANK
	REAL*8 WERF_EXP,IBOUND,IMIN
C
C PROGDESC is a variable use to confirm that the scratch block is not
C being used by some other routine.
C
	PROGDESC=1246719.0D0		!Must be unique
C
	IF(ND .GT. NV)THEN
	  WRITE(6,*)'Error in FG_HAM - NV not large enough'
	  STOP
	END IF
C
C Zero common block. There is currently 27 vectors in the common block.
C TA must be the first vector, and dCHIdR the last.
C
	dCHIdR(NV-1)=1.0D0
	dCHIdR(NV)=1.0D0
	I=(NV*27)-1
	CALL DP_ZERO(TA,I)
	IF(dCHIdR(NV-1) .NE. 0 .AND. dCHIdR(NV) .NE. 1)THEN
	  WRITE(6,*)'Error in zeroing SCRATCH block in FG_HAM'
	  STOP
	ELSE
	  dCHIdR(NV)=0.0D0
	END IF
C
C
C
C DNU is the frequency bandwidth over which we integrate to obtain
C JINT and is determined by the maximum expansion velocity of the
C atmosphere. We define DNU with out the factor FL.
C
	DNU=3.33564D-06*V(1)*2.0D0
C
C Save the line intensity in the large frequency band so that it can be
C used when we correct for electron scattering. If JNU has not been determined
C by a call to MOMJBAR, we set it equal to the continuum intensity.
C
	IF( JNU(1,NLF+1) .EQ. 0 .AND. JNU(ND,NLF+1) .EQ. 0 )THEN
	  DO I=1,ND
	    JPREV(I)=JCONT(I)*DNU
	  END DO
	ELSE
	  DO I=1,ND
	    JPREV(I)=JNU(I,NLF+1)
	  END DO
	END IF
C
C Zero intenisty matrices.
C
	J=(NLF+1)*ND
	CALL DP_ZERO(JNU,J)
	CALL DP_ZERO(HNU,J)
	CALL DP_ZERO(KNU,J)
	CALL DP_ZERO(NNU,J)
C
C Zero boundary condition vectors.
C
	DO ML=1,NLF+1
	  HBC(ML)=0.0D0
	  NBC(ML)=0.0D0
	  IN_HBC(ML)=0.0D0
	END DO
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
	      TCHI(I)=CHI(I)+CHIL(I)*MID_PRO
	      SOURCE(I)=(ETA(I)+ETAL(I)*MID_PRO)/TCHI(I)
	    END DO
	    CALL QKIM(Q,QH,GAM,GAMH,TCHI,PF,ML,NI,NLF)
	    IF(DIF .AND. LS .LE. NC)THEN
	      T1=CHI(ND)+CHIL(ND)*PROF(ML)
	      DBC=DBB*DSQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/TCHI(ND)
	1            *( 1.0D0+Q(ND)*TCHI(ND)*(1.0D0/T1-1.0D0/OLDCHI) )
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
C  IMIN=I(incident)=ETAL(1)/CHIL(1)*(1.0D0-WERF_EXP)+IBOUND*WERF_EXP
C
	    IF(THK_LINE .AND. ML .NE. 1)THEN
	      WERF_EXP=DEXP( 1.0D-15*CHIL(1)/FL/GAM(1)*
	1                      0.5D0*(WERFC(ML)+WERFC(ML-1)) )
	      XM(1)=(CHI(1)/TCHI(1)*WERF_EXP-1.0D0)*ETAL(1)/CHIL(1)
	1              -IBOUND*CHI(1)/TCHI(1)*WERF_EXP
	      IMIN=IBOUND*WERF_EXP + ETAL(1)/CHIL(1)*(1.0D0-WERF_EXP)
	    ELSE
	      XM(1)=-IBOUND
	      IMIN=IBOUND
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
C Solve for the radiation field along ray for this frequency.
C
	    CALL THOMAS(TA,TB,TC,AV,NI,1)
C
	    DO I=1,NI
	      IF(AV(I) .LT. 0)AV(I)=SOURCE(I)
	    END DO
C
C Update CV vector. NB Sign of GB is reveresd to that in CMFJBAR.
C
	    DO I=1,NI-1
	      CV(I)=GB(I)*(AV(I+1)-AV(I)+AVM1(I+1)-AVM1(I))+H(I)*CVM1(I)
	    END DO
C
C Update the Moments of the intensity.
C
	    DO I=1,NI
	      JNU(I,ML)=JNU(I,ML)+JQW(I,LS)*AV(I)
	      KNU(I,ML)=KNU(I,ML)+KQW(I,LS)*AV(I)
	    END DO
C
	    DO I=1,NI-1
	      HNU(I,ML)=HNU(I,ML)+HQW(I,LS)*CV(I)
	      NNU(I,ML)=NNU(I,ML)+NQW(I,LS)*CV(I)
	    END DO
C
C The quadrature weights used to compute HBC and NBC are crude.
C Their accuracy needs to be tested. Note that V=AV(1)-MIN.
C
	    HBC(ML)=HBC(ML)+JQW(1,LS)*(AV(1)-IMIN)*( Z(1)/R(1) )
	    NBC(ML)=NBC(ML)+JQW(1,LS)*(AV(1)-IMIN)*( (Z(1)/R(1)) )**3
	    IF(NI .EQ. ND)THEN
	      IN_HBC(ML)=IN_HBC(ML)+JQW(ND,LS)*( AV(ND)-(AV(ND)-AV(ND-1))
	1       /DTAU(ND-1) )*Z(ND)/R(ND)
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
	    IF(LINE_BL .AND. ML .EQ. 1)THEN
	      DO I=1,NI-1
	        CVCONT(I)=CV(I)
	        AVCONT(I)=AV(I)
	      END DO
	      AVCONT(NI)=AV(NI)
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
	      SOURCE(I)=( DNU*ETA(I)+ESEC(I)*(JPREV(I)
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
	    TA(1)=0.0D0
	    TC(1)=1./DTAU(1)
	    TB(1)=-1.0D0-TC(1)
	    DO I=2,NI-1
	      TA(I)=TC(I-1)
	      TC(I)=1.0D0/DTAU(I)
	      TB(I)=-0.5D0*(DTAU(I-1)+DTAU(I))-TA(I)-TC(I)
	      XM(I)=-0.5D0*( SOURCE(I)+Q(I) )*( DTAU(I-1)+DTAU(I) )
	1               -QH(I)+QH(I-1)
	    END DO
C
	    IF(LS .LE. NC .AND. DIF)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=TC(NI-1)
	      XM(NI)=DBC			!DNU include in definition of DBC.
	    ELSE IF(LS .GT. NC)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-TA(NI)+DTAU(NI-1)*0.5D0
	      XM(NI)=0.5D0*DTAU(NI-1)*(SOURCE(NI)+Q(NI))-QH(NI-1)
	    ELSE
	      TA(NI)=-TC(NI-1)
	      TB(NI)=1-TA(NI)
	      XM(NI)=IC*DNU+Q(NI)
	    END IF
	    TC(NI)=0.0
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
	      JNU(I,NLF+1)=JNU(I,NLF+1)+JQW(I,LS)*XM(I)
	      KNU(I,NLF+1)=KNU(I,NLF+1)+KQW(I,LS)*XM(I)
	    END DO
C
	    DO I=1,NI-1
	      HNU(I,NLF+1)=HNU(I,NLF+1)+HQW(I,LS)*CV(I)
	      NNU(I,NLF+1)=NNU(I,NLF+1)+NQW(I,LS)*CV(I)
	    END DO
C
C Note that the incident intensity is IBOUND*DNU.
C
	    HBC(NLF+1)=HBC(NLF+1)+JQW(1,LS)*(XM(1)-IBOUND*DNU)*Z(1)/R(1)
	    IF(NI .EQ. ND)THEN
	      IN_HBC(NLF+1)=IN_HBC(NLF+1)+JQW(ND,LS)*
	1     (XM(ND)-(XM(ND)-XM(ND-1))/DTAU(ND-1))*Z(ND)/R(ND)
	    END IF
C
	  END IF		!End if(LINE_BL)
C
	END DO			!End do LS
C
C
C Compute the Eddington F and G factors, which are defined by
C F=K/J and G=N/H. The F and G factors are returned in KNU and NNU
C respectively.
C
	J=NLF
	IF(LINE_BL)J=NLF+1
	DO ML=1,J
	  DO I=1,ND-1
	    KNU(I,ML)=KNU(I,ML)/JNU(I,ML)
	    NNU(I,ML)=NNU(I,ML)/HNU(I,ML)
	  END DO
	  KNU(ND,ML)=KNU(ND,ML)/JNU(ND,ML)
	END DO
C
C Compute the boundary Eddington factors.
C
	DO ML=1,J
	  HBC(ML)=HBC(ML)/JNU(1,ML)
	  NBC(ML)=NBC(ML)/JNU(1,ML)
	  IF(ML .EQ. NLF+1)THEN
	    IN_HBC(ML)=IN_HBC(ML)/(2.0D0*JNU(ND,ML)-IC*DNU)
	  ELSE
	    IN_HBC(ML)=IN_HBC(ML)/(2.0D0*JNU(ND,ML)-IC)
	  END IF
	END DO
C
C
	IF(LINE_BL)THEN
C
C Evaluate HBLANK. Scale HBLANK to allow for the fact that the LFQW
C dont necessarily add to v*(PF(1)-PF(NLF)) (since LFQW is normalized so
C that integral over the line profile is unity).
C
	  HBLANK=0.0D0
	  SCALE=0.0D0
	  DO ML=1,NLF
	    HBLANK=HBLANK+LFQW(ML)*HNU(1,ML)
	    SCALE=SCALE+LFQW(ML)
	  END DO
	  SCALE=1.0D+15*FL*(PF(1)-PF(NLF))/SCALE
	  HBLANK=HBLANK*SCALE
C
C Evaluate HBLANK which is deined by Int{Jv dv} over the WHOLE line.
C Note that we need to keep JINT since this is used in the electron
C scattering integral. The frequency factor of 10^15 is included in the
C HBLANK definition.
C
	  T1=1.0D+15*FL
	  HBLANK=HBLANK+HNU(1,NLF+1)*T1
C
C Evaluate the line EW. The units are Angstroms. Also evaluate
C the continuum intensity ( Jys/kpc/kpc ). Note that H is
C defined midway between R(1) and R(2).
C
	  T1=( (PF(1)-PF(NLF))+DNU )*FL*1.0D+15
	  EW=2.99794D-12*( HBLANK-HNU(1,1)*T1 )/HNU(1,1)/FL/FL
	  CONT_INT=13.19868*HNU(1,1)*( (R(1)+R(2))**2 )/4.0D0
	END IF
C
C
C
	IF(PROGDESC .NE. 1246719.0D0)THEN
	  WRITE(6,*)'Error - Scratch Block Corrupted in FG_HAM'
	  STOP
	END IF
C
	RETURN
	END
