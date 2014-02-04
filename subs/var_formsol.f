C
C This subroutine computes the variation of Jmn as a function of both
C the line and continuum emissivities and opacities. The solution is
C done in the (p,z) plane using linear differencing in frequency.
C
C Zero FQAF matrix and FQAFD vector. FQAF is of dimension (ND,ND,NM)
C and is equaly initially to dJ/dx where x=(chil,etal,chi,eta) respectively.
C FQAFD is initially equal to dJ/d(dT/dR). Note that dT/dR is only a function
C of the populations at the inner boundary.
C
	SUBROUTINE VAR_FORMSOL(ETA,CHI,ESEC,CHIL,ETAL,V,SIGMA,R,P,
	1                      FQAF,FQAFD,TX,TVX,KI,JQW,
	1                      PF,PROF,LFQW,WERFC,FL,DIF,DBB,IC,AMASS,
	1                      THK_LINE,THK_CONT,
	1                      NLF,NC,NP,ND,NM,METHOD)
	IMPLICIT NONE
C
C Altered 28-Oct-1996 - Bug fix: COS converted back to ACOS in TOR expression.
C Altered 28-May-1996 - Scratch block removed.
C                       Dynamic allocation for work vectors.
C                       Call to DP_ZERO removed.
C                       IONE inserted in call to SIMPTH
C                       Genrical calls for EXP.
C Created 12-May-1989 - Based on part of [JDH.FINAL]LINEGEN.FOR and FORMSOL.
C
	INTEGER NLF,NC,NP,ND,NM
	REAL*8 ETA(ND),CHI(ND),ESEC(ND),CHIL(ND),ETAL(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND),P(NP)
	REAL*8 FQAF(ND,ND,NM),FQAFD(ND)
	REAL*8 TX(ND,ND,NM),TVX(ND-1,ND,NM),KI(ND,3,NM)
	REAL*8 JQW(ND,NP)
	REAL*8 PF(NLF),PROF(NLF),LFQW(NLF),WERFC(NLF)
	REAL*8 DBB,IC,AMASS,FL
	LOGICAL DIF,THK_CONT,THK_LINE
	CHARACTER*(*) METHOD			!Not yet implemented
C
C Wrok vectors.
C
	REAL*8 TA(ND),TB(ND),TC(ND),AV(ND),CV(ND),DTAU(ND),Z(ND)
	REAL*8 TCHI(ND),XM(ND),SOURCE(ND),U(ND),VB(ND),VC(ND)
	REAL*8 GB(ND),H(ND),GAM(ND),GAMH(ND),Q(ND),QH(ND)
	REAL*8 RKB(ND),RKC(ND),BAR_WTS(ND),DEPTH_PRO(ND)
	REAL*8 AVM1(ND),CVM1(ND),TXD(ND),TVXD(ND)
C
C Local variables.
C
	INTEGER, PARAMETER :: IONE=1
C
	LOGICAL MLNE1
	INTEGER I,LS,ML,NI
	REAL*8 OLDCHI,T1,T2,DBC,DBC_ON_DBB,TOR,IBOUND,WERF_EXP
C
	FQAF(:,:,:)=0.0D0               !NM,ND,ND
	FQAFD(:)=0.0D0         		!ND
C
C 
C
C Enter impact parameter loop.
C
	DO LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)NI=ND
	  CALL ZALONGP(R,Z,P(LS),NI)
	  CALL GAMMA(GAM,GAMH,SIGMA,Z,R,V,ND,NI)
C
C Determine boundary condition for continuum intensity.
C
	  IF(THK_CONT)THEN
	    IF(P(LS) .GT. 0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796-ACOS(P(LS)/R(1)))/P(LS)
	    ELSE
	      TOR=CHI(1)*R(1)
	    END IF
	    IBOUND=ETA(1)*(1.0D0-EXP(-TOR))/CHI(1)
	  ELSE
	    TOR=0.0D0
	    IBOUND=0.0D0
	  END IF
C
C Zero AV and CV vectors.
C
	  AV(1:NI)=0.0D0
	  CV(1:NI)=0.0D0
C
C set MLNE1 to false for first frequency.
	  MLNE1=.FALSE.
C
C 
C
C Perform integration for each frequency in turn.
C
	  OLDCHI=CHI(NI)
	  DO ML=1,NLF
C
C Store previous U,V intensity values (also zero's AVM1 and CVM1
C for first freqency)
C
	    DO I=1,NI
	      AVM1(I)=AV(I)
	      CVM1(I)=CV(I)
	    END DO
C
	    DO I=1,NI
	      DEPTH_PRO(I)=PROF(ML)
	      TCHI(I)=CHI(I)+PROF(ML)*CHIL(I)
	      SOURCE(I)=(ETA(I)+ETAL(I)*PROF(ML))/TCHI(I)
	    END DO
	    CALL QKIM(Q,QH,GAM,GAMH,TCHI,PF,ML,NI,NLF)
	    IF(DIF .AND. LS .LE. NC)THEN
	      DBC_ON_DBB=DSQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/TCHI(ND)
	1                   *(1.0D0+Q(NI)*(1.0D0-TCHI(NI)/OLDCHI))
	      DBC=DBB*DBC_ON_DBB
	    END IF
	    CALL TAU(DTAU,TCHI,Z,NI)
	    CALL TUVGHD(TA,TB,TC,U,VB,VC,GB,H,Q,QH,DTAU,DIF,LS,NC,NI)
	    CALL XVECD(DTAU,SOURCE,XM,DIF,DBC,IC,LS,NC,ND,NI)
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
C Solve for the radiation field along ray for this frequency.
C
	    CALL THOMAS(TA,TB,TC,AV,NI,1)
C
C Update Flux vector.
C
	    DO I=1,NI-1
	      CV(I)=GB(I)*(AV(I)-AV(I+1))+H(I)*CV(I)
	    END DO
C
C 
C
C The u and v components of the radiation field have been found. We can thus
C begin the variation computation.
C
C Compute d{non-radiation field}/dchi matrix.
C
	    T1=1.0D-15/FL/GAM(1)*WERFC(ML)
	    CALL KIVARNM(KI,RKB,RKC,DEPTH_PRO,DTAU,Z,Q,QH,
	1                AV,AVM1,CVM1,TCHI,
	1                SOURCE,ETAL,CHIL,CHI,T1,LS,NC,NI,NM)
C
C Include continuum term in outer boundary condition.
C
	    IF(THK_LINE)THEN
	      KI(1,1,1)=KI(1,1,1) + IBOUND*CHI(1)*WERF_EXP *
	1                            DEPTH_PRO(I)/TCHI(1)/TCHI(1)
	      IF(NM .EQ. 4)KI(1,1,3)=KI(1,1,3) - IBOUND*WERF_EXP/TCHI(1)
	    END IF
C
C Treat diffusion boundary conndition.
C
	    IF(DIF .AND. LS .LE. NC)THEN
	      T1=0
	      IF(MLNE1)T1=PROF(ML-1)
	      T2=0.5D0*(Z(ND-1)-Z(ND))*(AV(ND)-AV(ND-1))/DTAU(ND-1)/DTAU(ND-1)
C
C We note that KI(ND,1,1) and KI(ND,1,3) are unaltered.
C
	      IF(NM .EQ. 4)THEN
	        KI(ND,2,3)=T2+DBC*( (Q(ND)*TCHI(ND)-1.0D0)/TCHI(ND)
	1        +Q(ND)/(OLDCHI**2) )
	      END IF
	      KI(ND,2,1)=T2*DEPTH_PRO(ND) +
	1         DBC*( (Q(ND)*TCHI(ND)-1.0D0)/TCHI(ND) *
	1         DEPTH_PRO(ND)+Q(ND)*T1/(OLDCHI**2) )
	    END IF
C
C Evaluat the intensity variations.
C                                  TX=d(u).d(chil,etal,chi,eta)
C                          and     TVX=d(v).d(chil,etal,chi,eta)
C
	    CALL UPDATEU(TX,TVX,KI,U,VB,VC,TA,TB,TC,MLNE1,NI,NM)
	    MLNE1=.TRUE.
	    CALL UPVNOT(TVX,TX,RKB,RKC,DEPTH_PRO,GB,H,NI,NM)
C
C Compute the "weights" to increment d{mean intensity}d{ , , , } arrays.
C
	    DO I=1,NI
	      BAR_WTS(I)=JQW(I,LS)*LFQW(ML)*DEPTH_PRO(I)
	    END DO
C
C Evaluate Diffusion variation.
C
	    IF(DIF .AND. LS .LE. NC)THEN
	      IF(ML .EQ. 1)THEN
	        TXD(1:NI)=0.0D0
	        TVXD(1:NI)=0.0D0
	      END IF
	      DO I=2,NI-1
	        TXD(I)=U(I)*TXD(I)-VB(I)*TVXD(I-1)-VC(I)*TVXD(I)
	      END DO
	      TXD(1)=U(1)*TXD(1)
	      TXD(NI)=DBC_ON_DBB+U(NI)*TXD(NI)-VB(NI)*TVXD(NI-1)
C
	      CALL SIMPTH(TA,TB,TC,TXD,NI,IONE)
C
	      DO I=1,NI-1
	        TVXD(I)=GB(I)*(TXD(I)-TXD(I+1))+H(I)*TVXD(I)
	        FQAFD(I)=FQAFD(I)+BAR_WTS(I)*TXD(I)
	      END DO
	      FQAFD(ND)=FQAFD(ND)+BAR_WTS(ND)*TXD(ND)
	    END IF	    	    !DIFF END
C
	    CALL MULT2D(FQAF,BAR_WTS,TX,ND,NI,NM)
	    OLDCHI=TCHI(NI)
	  END DO
	END DO
C
	RETURN
	END
