C
C Routine to compute the Eddington F and G factors for line transfer
C in the comoving-frame. Based on CMFJBAR. First order differencing
C in frequency is used.
C
	SUBROUTINE FG_COMP(ETA,CHI,ESEC,JCONT,
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
C Altered 28-Oct-1996 - Bug fix: COS corrected back to ACOS in TOR expression.
C Altered 24-May-1996 - SCRTEMP common block removed.
C                       CALLs to DP_ZERO removed.
C                       No longer a limit on ND.
C                       IONE used in CALL to THOMAS
C                       Generical calls for EXP, COS.
C
C Altered 15-Jun-1989 - Outer boundary condition changed. HBC and NBC must now
C                       be dimensioned [ 3, (NLF+1) ]. Altered so that SL
C                       appears explicitliy in bc. Boundary condition for
C                       line blanketing section unchaanged.
C
C Altered 8-Jun-1989 -  Continuum boundary corrected (wasn't dividing by CHI).
C                       Line boundary condition improved. Now compute dI-/dnu
C                       numerically at outer boundary. Its computation
C                       analytically introdueced an artificial source term into
C                       the EW computation.
C
C Altered 12-May-1989 - Cleaning. Thick boundary condition fixed in EW
C                       section.
C Altered  8-May-1989 - (and April). Cleaning and bug fixes.
C                       EW computation installed.
C Created 11-JAN-1989
C
	INTEGER NLF,NC,NP,ND
	REAL*8 ETA(ND),CHI(ND),ESEC(ND),JCONT(ND),CHIL(ND),ETAL(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND),P(NP)
C
	REAL*8 JQW(ND,NP),HQW(ND,NP),KQW(ND,NP),NQW(ND,NP)
	REAL*8 JNU(ND,NLF+1),HNU(ND,NLF+1),KNU(ND,NLF+1),NNU(ND,NLF+1)
	REAL*8 IN_HBC(NLF+1),HBC(3,NLF+1),NBC(3,NLF+1)
C
	REAL*8 PF(NLF),PROF(NLF),LFQW(NLF),WERFC(NLF)
	REAL*8 DBB,IC,FL,EW,CONT_INT
	CHARACTER*6 METHOD
	LOGICAL DIF,LINE_BL,THK_CONT,THK_LINE
C
	REAL*8 TA(ND),TB(ND),TC(ND),AV(ND),CV(ND),DTAU(ND),Z(ND)
	REAL*8 TCHI(ND),XM(ND),SOURCE(ND),U(ND),VB(ND),VC(ND)
	REAL*8 GB(ND),H(ND),GAM(ND),GAMH(ND),Q(ND),QH(ND)
	REAL*8 JPREV(ND),AVCONT(ND),CVCONT(ND),dCHIdR(ND)
C
C Local variables.
C
	INTEGER, PARAMETER :: IONE=1
	INTEGER I,J,LS,ML,NI,NIEXT
	REAL*8 OLDCHI,DBC,SCALE,HBLANK,IBOUND,IMIN,IMINOLD,IMIN_dSL
	REAL*8 T1,TOR,DNU,WERF_EXP
C
	TA(:)=0.0D0;      TB(:)=0.0D0;     TC(:)=0.0D0;  AV(:)=0.0D0;
        CV(:)=0.0D0;      DTAU(:)=0.0D0;   Z(:)=0.0D0
	TCHI(:)=0.0D0;    XM(:)=0.0D0;     SOURCE(:)=0.0D0
	U(:)=0.0D0;       VB(:)=0.0D0;     VC(:)=0.0D0
	GB(:)=0.0D0;      H(:)=0.0D0;      GAM(:)=0.0D0; GAMH(:)=0.0D0
	Q(:)=0.0D0;       QH(:)=0.0D0;     JPREV(:)=0.0D0
	AVCONT(:)=0.0D0;  CVCONT(:)=0.0D0; dCHIdR(:)=0.0D0
 
C
C 
C
C DNU is the frequency bandwidth over which we integrate to obtain
C JINT and is determined by the maximum expansion velocity of the
C atmosphere. We define DNU with out the factor FL. Must be identical
C to that in MOMJBAR.
C
	DNU=3.33564D-06*V(1)*2.0D0
C
C Save the line intensity in the large frequency band so that it can be
C used when we correct for electron scattering. If JNU has not been determined
C by a call to MOMJBAR, we set it equal to the continuum intensity.
C
	IF( JNU(1,NLF+1) .EQ. 0.0D0 .AND. JNU(ND,NLF+1) .EQ. 0.0D0 )THEN
	  DO I=1,ND
	    JPREV(I)=JCONT(I)*DNU
	  END DO
	ELSE
	  DO I=1,ND
	    JPREV(I)=JNU(I,NLF+1)
	  END DO
	END IF
C
C Zero intenisty matrices. All are dimensions (ND,NLF+1)
C
	JNU(:,:)=0.0D0
	HNU(:,:)=0.0D0
	KNU(:,:)=0.0D0
	NNU(:,:)=0.0D0
C
C Zero boundary condition vectors. Dimensioned (3,NLF+1)
C
	HBC(:,:)=0.0D0
	NBC(:,:)=0.0D0
C
	DO ML=1,NLF+1
	  IN_HBC(ML)=0.0D0
	END DO
C
C Enter loop to perform integration along each ray.
C
	DO LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)NI=ND
	  IF(NI .EQ. ND)THEN
	    NIEXT=ND
	  ELSE
	    NIEXT=NI+1
	  END IF
C
C Zero AV and CV vectors.
C
	  AV(1:NI)=0.0D0
	  CV(1:NI)=0.0D0
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
C
C We evaluate up to NIEXT for DERIVCHI. By setting PF(1)=0 when
C evaluating TCHI and SOURCE we ensure a pure continuum calculation
C for the first frequency.
C
	    IF(ML .EQ. 1)THEN
	      T1=0.0D0
	    ELSE
	      T1=PROF(ML)
	    END IF
	    DO I=1,NIEXT
	      TCHI(I)=CHI(I)+CHIL(I)*T1
	      SOURCE(I)=(ETA(I)+ETAL(I)*T1)/TCHI(I)
	    END DO
	    CALL QKIM(Q,QH,GAM,GAMH,TCHI,PF,ML,NI,NLF)
	    IF(DIF .AND. LS .LE. NC)THEN
	      DBC=DBB*DSQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/TCHI(ND)
	1   *(1.0D0+Q(NI)*(1.0D0-TCHI(NI)/OLDCHI))
	    END IF
C
	    IF(METHOD .EQ. 'ZERO')THEN
	      CALL TAU(DTAU,TCHI,Z,NI)
	    ELSE
	      CALL DERIVCHI(dCHIdR,TCHI,R,NIEXT,METHOD)
	      CALL NORDTAU(DTAU,TCHI,Z,R,dCHIdR,ND)
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
C If THK_LINE is TRUE, we adopt the SOBOLEV approximation for
C the incident intensity at the outer boundary. The incident intensity
C is both angle and frequency independent. The reason for the ratio
C CHI/TCHI is related to the the first order equation for U.
C
C Note that WERFC=-0.5*ERFC where ERFC is the complementary error
C function. We use WERF_EXP because the value for ML=NLF is
C required in the line blanketing section.
C
C	    I(incident)=ETAL(1)/CHIL(1)*(1.0D0-WERF_EXP)+IBOUND*WERF_EXP
C
C Now compute dI-/dnu numerically at outer boundary. Its computation
C analytically introduced an artificial source term into the EW computation.
C Could symetrize the line boundary condition by replacing the WERF_EXP
C calculation by the following statements.
C
C	    WERF_EXP=WERFC(ML)
C	    IF(WERFC(ML) .LT. -0.5D0)WERF_EXP=-1.0D0-WERFC(ML)
C	    WERF_EXP=EXP( 1.0D-15*CHIL(1)/FL/GAM(1)*WERF_EXP )
C
	    IF(THK_LINE)THEN
	      WERF_EXP=EXP( 1.0D-15*CHIL(1)/FL/GAM(1)*WERFC(ML) )
	      IF(ML .EQ. 1)THEN
	        IMINOLD=0.0D0
	      ELSE
	        IMINOLD=IMIN
	      END IF
	      IMIN=ETAL(1)/CHIL(1)*(1.0D0-WERF_EXP)+IBOUND*WERF_EXP
	      IMIN_dSL=(1.0D0-WERF_EXP)
	      XM(1)=Q(1)*(IMINOLD-IMIN)-IMIN
	    ELSE
	      XM(1)=-IBOUND
	      IMIN=IBOUND
	      IMIN_dSL=0.0D0
	    END IF
C
C Update AV matrix.
C
	    AV(1)=XM(1)+U(1)*AV(1)
	    DO I=2,NI-1
	      AV(I)=XM(I)+( U(I)*AV(I)-(VB(I)*CV(I-1)+VC(I)*CV(I)) )
	    END DO
	    AV(NI)=XM(NI)+( U(NI)*AV(NI)-VB(NI)*CV(NI-1) )
C
C Solve for the radiation field along ray for this frequency.
C
	    CALL THOMAS(TA,TB,TC,AV,NI,IONE)
C
C Update C vector (i.e. flux variable).
C
	    DO I=1,NI-1
	      CV(I)=GB(I)*(AV(I)-AV(I+1))+H(I)*CV(I)
	    END DO
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
C Their accuracy needs to be tested. Note that V=AV(1)-IMIN.
C
	    T1=Z(1)/R(1)
	    HBC(1,ML)=HBC(1,ML)+JQW(1,LS)*AV(1)*T1
	    HBC(2,ML)=HBC(2,ML)+JQW(1,LS)*IMIN*T1
	    HBC(3,ML)=HBC(3,ML)+JQW(1,LS)*IMIN_dSL*T1
C
	    T1=T1**3
	    NBC(1,ML)=NBC(1,ML)+JQW(1,LS)*AV(1)*T1
	    NBC(2,ML)=NBC(2,ML)+JQW(1,LS)*IMIN*T1
	    NBC(3,ML)=NBC(3,ML)+JQW(1,LS)*IMIN_dSL*T1
C
	    IF(NI .EQ. ND)THEN
	      IN_HBC(ML)=IN_HBC(ML)+
	1       JQW(ND,LS)*( AV(ND)-(AV(ND)-AV(ND-1))
	1       /DTAU(ND-1) )*Z(ND)/R(ND)
	    END IF
	    OLDCHI=TCHI(NI)
C
C Store the continuum mean and flux intensities along this ray. They
C are required for the line blanketing calculation.
C
	    IF(ML .EQ. 1)THEN
	      DO I=1,NI
	        AVCONT(I)=AV(I)
	        CVCONT(I)=CV(I)
	      END DO
	    END IF
	  END DO		!end do ML
C 
C
C This section of the routine provides a means of estimating the
C line blanketing for selected lines. Once where on the red side of the
C line profile the opacity is constant, hence we can integrate
C the transfer equation with respect to frequency. This leaves
C us with a simple "continuum like" transfer problem. To
C do the frequency integration it is also necessary to assume that
C the source function is frequency independent. We have thus assumed
C that the line photons are scattered coherently. This should be an
C excellent approximation since the continuum transfer of the line
C photons should not be strongly influenced by where the photons are.
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
C Q and QH are temporary vectors which arrises from the frequency integration
C of the d_nu term in the transfer equation. They are defined differently
C in the line section of the program.
C
	    DO I=1,NI-1
	      Q(I)=GAM(I)*( AV(I)-AVCONT(I) )/CHI(I)
	      QH(I)=2.0D0*GAMH(I)*( CV(I)-CVCONT(I) )/(CHI(I)+CHI(I+1))
	    END DO
	    Q(NI)=GAM(NI)*( AV(NI)-AVCONT(NI) )/CHI(NI)
C
C For the line, it is very difficult to estimate the incident
C line intensity once where outside the core. We thus assume that
C it is given by the continuum intensity = IBOUND. NB - The term
C added for the case of a thick line is GAM''*(I(blue) - Icont).
C
	    XM(1)=-Q(1)-IBOUND*DNU
	    IF(THK_LINE)XM(1)=XM(1)+GAM(1)*
	1            (ETAL(1)/CHIL(1)-IBOUND)*(1.0D0-WERF_EXP)/CHI(1)
	    TA(1)=0.0D0
	    TC(1)=1.0D0/DTAU(1)
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
	      TB(NI)=1.0D0-TA(NI)
	      XM(NI)=IC*DNU+Q(NI)
	    END IF
	    TC(NI)=0.0D0
C
C Solve the tridiagonal system of equations.
C
	    CALL THOMAS(TA,TB,TC,XM,NI,IONE)
C
C Compute the flux
C
	    DO I=1,NI-1
	      CV(I)=(XM(I+1)-XM(I))/DTAU(I)+QH(I)
	    END DO
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
	    HBC(1,NLF+1)=HBC(1,NLF+1)+
	1                 JQW(1,LS)*(XM(1)-IBOUND*DNU)*Z(1)/R(1)
	    IF(NI .EQ. ND)THEN
	      IN_HBC(NLF+1)=IN_HBC(NLF+1)+JQW(ND,LS)*
	1     (XM(ND)-(XM(ND)-XM(ND-1))/DTAU(ND-1))*Z(ND)/R(ND)
	    END IF
	  END IF		!End if(LINE_BL)
C 
	END DO			!End do LS
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
	  HBC(1,ML)=HBC(1,ML)/JNU(1,ML)
	  NBC(1,ML)=NBC(1,ML)/JNU(1,ML)
	  IF(ML .EQ. NLF+1)THEN
	    IN_HBC(ML)=IN_HBC(ML)/(2.0D0*JNU(ND,ML)-IC*DNU)
	  ELSE
	    IN_HBC(ML)=IN_HBC(ML)/(2.0D0*JNU(ND,ML)-IC)
	  END IF
	END DO
C
	IF(LINE_BL)THEN
C
	  HBLANK=0.0D0
	  DO ML=1,NLF
	    HBLANK=HBLANK+LFQW(ML)*HNU(1,ML)
	  END DO
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
C Evaluate JBLANK which is deined by Int{Jv dv} over the WHOLE line.
C HBLANK is similarly defined. WE first scale JBLANK and HBLANK for
C the reasons givem above. The frequency factor of
C 10^15 is included in the JBLANK (and HBLANK) definition.
C
	  T1=1.0D+15*FL
	  HBLANK=HBLANK*SCALE+HNU(1,NLF+1)*T1
C
C Evaluate the line EW. The units are Angstroms. Also evaluate
C the continuum intensity ( Jys/kpc/kpc ). Note that H is
C defined midway between R(1) and R(2).
C
	  T1=( (PF(1)-PF(NLF))+DNU )*FL*1.0D+15
	  EW=2.99794D-12*( HBLANK-HNU(1,1)*T1 )/HNU(1,1)/FL/FL
	  CONT_INT=13.19868D0*HNU(1,1)*( (R(1)+R(2))**2 )/4.0D0
	ELSE
	  EW=0.0D0
	  CONT_INT=0.0D0
	END IF
C
	RETURN
	END
