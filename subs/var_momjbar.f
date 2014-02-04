C
C This subroutine computes the variation of Jmn as a function of both
C the line and continuum emissivities and opacities. The solution is
C done in the (p,z) plane using linear differencing in frequency.
C
C Zero FQAF matrix and FQAFD vector. FQAF is of dimension (ND,ND,NM)
C and is equally initially to dJ/dx where x=(chil,etal,chi,eta) respectively.
C FQAFD is initially equal to dJ/d(dT/dR). Note that dT/dR is only a function
C of the populations at the inner boundary.
C
C The particular choice of the outer boundary condition adopted is irrelevant
C for this routine. Such information is incorporated by the outer boundary
C Eddington factors HBC, and NBC.
C
C The program computes r^2 JNU and r^2 HNU.
C
	SUBROUTINE VAR_MOMJBAR(ETA,CHI,ESEC,CHIL,ETAL,
	1                  V,SIGMA,R,
	1                  FQAF,FQAFD,TX,TVX,
	1                  KI,WORKMAT,RHS_dHdCHI,
	1                  F,G,HBC,IN_HBC,NBC,
	1                  PF,PROF,LFQW,FL,DIF,DBB,IC,METHOD,
	1                  NLF,ND,NM)
	IMPLICIT NONE
C
C Altered 02-Jul-1998 : LUER and ERROR_LU installed.
C Altered 22-Sep-1994 : After EDDLINE_VAR, KI is updated by T1 to allow for
C                       the boundary condition. For ML=1, the computation of T1
C                       accessed values outside the bounds of HBC and NBC.
C                       Had no effect, since accessed terms multiplied by zero.
C Altered 15-Jun-1989 - Boundary condition changed, HBC and NBC are now
C                       matrices dimensioned [3,NLF+1]. EDDLINE_VAR
C                       not altered.
C Altered 01-Jun-1989 - THETA, NC and NP removed from call.
C
C Created 16-May-1989 - Based on part of [JDH.EDDHAM]MOMHAM.FOR
C                                    and [JDH.RAYLINE]VAR_FORMSOL.FOR
C
	INTEGER NLF,ND,NM
	REAL*8 ETA(ND),CHI(ND),ESEC(ND)
	REAL*8 CHIL(ND),ETAL(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND)
C
C Varaiation arrays and vectors.
C
	REAL*8 FQAF(ND,ND,NM),FQAFD(ND)
	REAL*8 TX(ND,ND,NM),TVX(ND-1,ND,NM)
	REAL*8 KI(ND,ND,NM),WORKMAT(ND,ND),RHS_dHdCHI(ND-1,ND)
C
C Radiation field variables.
C
	REAL*8 F(ND,NLF+1),G(ND,NLF+1)
	REAL*8 HBC(3,NLF+1),NBC(3,NLF+1),IN_HBC(NLF+1)
C
C Profile information
C
	REAL*8 PF(NLF),PROF(NLF),LFQW(NLF)
C
	REAL*8 DBB,IC,FL
	CHARACTER*6 METHOD
	LOGICAL DIF
C
	INTEGER NV
	PARAMETER (NV=100)
	COMMON /SCRATCH/  PROGDESC,TA,TB,TC,DTAU,TCHI,DTAUONQ,
	1                   XM,SOURCE,MIDF,Q,
	1                   JNU,HNU,JNUM1,HNUM1,TXD,TVXD,
	1                   BAR_WTS,DEPTH_PRO,
	1                   VB,VC,HU,HL,HS,
	1                   GAM,GAMH,W,WPREV,PSI,PSIPREV
C
	REAL*8 TA(NV),TB(NV),TC(NV),DTAU(NV),TCHI(NV),DTAUONQ(NV)
	REAL*8 XM(NV),SOURCE(NV),MIDF(NV),Q(NV)
	REAL*8 JNU(NV),HNU(NV),JNUM1(NV),HNUM1(NV),TXD(NV),TVXD(NV)
	REAL*8 BAR_WTS(NV),DEPTH_PRO(NV)
	REAL*8 VB(NV),VC(NV),HU(NV),HL(NV),HS(NV)
	REAL*8 GAM(NV),GAMH(NV),W(NV),WPREV(NV)
	REAL*8 PSI(NV),PSIPREV(NV)
C
	REAL*8 PROGDESC	
	REAL*8, PARAMETER :: PROG_ID=2.22814021D+08 !Must be unique (VAR_MJBA)
C
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
C
C Local variables.
C
	LOGICAL ML_NE_ONE
	INTEGER I,J,K,ML
	REAL*8 T1,DNU,SRCEBND
C
C PROGDESC is a variable use to confirm that the scratch block is not
C being used by some other routine.
C
	PROGDESC=PROG_ID
	IF(ND .GT. NV)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in VAR_MOMJBAR - NV smaller than ND'
	  WRITE(LUER,*)'ND=',ND,'NV',NV
	  STOP
	END IF
C
C Zero common block. There are currently 29 vectors in the common block.
C TA must be the first vector, and PSIPREV the last.
C
	PSIPREV(NV-1)=1.0D0
	PSIPREV(NV)=1.0D0
	I=(NV*29)-1
	CALL DP_ZERO(TA,I)
	IF(PSIPREV(NV-1) .NE. 0 .OR. PSIPREV(NV) .NE. 1)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in zeroing SCRATCH block in VAR_MOMJBAR'
	  STOP
	ELSE
	  PSIPREV(NV)=0.0D0
	END IF
C
C 
C
C Zero relevant vectors and matrices.
C
	DO I=1,ND
	  JNU(I)=0.0D0
	  HNU(I)=0.0D0
	  JNUM1(I)=0.0D0
	  HNUM1(I)=0.0D0
	  TXD(I)=0.0D0
	  TVXD(I)=0.0D0
	END DO
C
	CALL DP_ZERO(FQAF, ND*ND*NM )
	CALL DP_ZERO(FQAFD, ND )
	CALL DP_ZERO(TX, ND*ND*NM )
	CALL DP_ZERO(TVX, (ND-1)*ND*NM )
C
C*****************************************************************************
C
	ML_NE_ONE=.FALSE.
	DO ML=1,NLF
C
C Compute the total opacity. Store opacity from previous frequency.
C Compute the total source function.
C
C                 *************************
C                 *************************
C It is currently assumed that ETA is the total continuum emissivity
C and hence already contains the continuum scattering term. This may
C need to be altered.
C                 *************************
C                 *************************
C
	  IF(ML .EQ. 1)THEN
	    DO I=1,ND
	      DEPTH_PRO(I)=0.0D0
	      TCHI(I)=CHI(I)
	      SOURCE(I)=ETA(I)/TCHI(I)
	    END DO
	  ELSE
	    DO I=1,ND
	      DEPTH_PRO(I)=PROF(ML)
	      TCHI(I)=CHI(I)+CHIL(I)*PROF(ML)
	      SOURCE(I)=(ETA(I)+ETAL(I)*PROF(ML))/TCHI(I)
	    END DO
	  END IF
C
C Compute the Q factors from F. Then compute optical depth scale.
C
	  CALL QFROMF(F(1,ML),Q,R,TA,TB,ND)	!TA work vector
	  DO I=1,ND
	    TA(I)=TCHI(I)*Q(I)
	  END DO
	  CALL DERIVCHI(TB,TA,R,ND,METHOD)
C
C We need to call d_DERIVCHI_dCHI to set the TRAP derivatives.
C
	  CALL d_DERIVCHI_dCHI(TB,TA,R,ND,METHOD)
	  CALL NORDTAU(DTAU,TA,R,R,TB,ND)
C
C 
C
C Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
C 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
C Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
C
	  IF(ML .NE. 1)THEN
	    DNU=PF(ML-1)-PF(ML)
C
C NB - By definition, G is defined at the mesh midpoints.
C
	    DO I=1,ND-1
	      GAMH(I)=2.0D0*3.33564D-06*(V(I)+V(I+1))/(R(I)+R(I+1))
	1         /DNU/( TCHI(I)+TCHI(I+1) )
	      W(I)=GAMH(I)*( 1.0D0+
	1                 0.5D0*(SIGMA(I)+SIGMA(I+1))*G(I,ML) )
	      WPREV(I)=GAMH(I)*( 1.0D0+
	1                 0.5D0*(SIGMA(I)+SIGMA(I+1))*G(I,ML-1) )
	    END DO
C
	    DO I=1,ND
	      GAM(I)=3.33564D-06*V(I)/R(I)/TCHI(I)/DNU
	    END DO
C
C PSIPREV is equivalent to the U vector of FORMSOL.
C
	    PSI(1)=GAM(1)*(HBC(1,ML)+NBC(1,ML)*SIGMA(1))
	    PSIPREV(1)=GAM(1)*(HBC(1,ML-1)+NBC(1,ML-1)*SIGMA(1))
	    SRCEBND=GAM(1)*( HBC(2,ML)+NBC(2,ML)*SIGMA(1)
	1                   -HBC(2,ML-1)-NBC(2,ML-1)*SIGMA(1) )
	    DO I=2,ND
	      DTAUONQ(I)=0.5D0*(DTAU(I)+DTAU(I-1))/Q(I)
	      PSI(I)=DTAUONQ(I)*GAM(I)*( 1.0D0+SIGMA(I)*F(I,ML) )
	      PSIPREV(I)=DTAUONQ(I)*GAM(I)*
	1                   ( 1.0D0+SIGMA(I)*F(I,ML-1) )
	    END DO
C
	  ELSE
	    DO I=2,ND
	      DTAUONQ(I)=0.5D0*(DTAU(I)+DTAU(I-1))/Q(I)
	    END DO
	    DO I=1,ND
	      GAMH(I)=0.0D0
	      W(I)=0.0D0
	      WPREV(I)=0.0D0
	      GAM(I)=0.0D0
	      PSI(I)=0.0D0
	      PSIPREV(I)=0.0D0
	    END DO
	    SRCEBND=0.0D0
	  END IF
C 
C
C Compute vectors used to compute the flux vector H.
C
	  DO I=1,ND-1
	    HU(I)=F(I+1,ML)*Q(I+1)/(1.0D0+W(I))/DTAU(I)
	    HL(I)=F(I,ML)*Q(I)/(1.0D0+W(I))/DTAU(I)
	    HS(I)=WPREV(I)/(1.0D0+W(I))
	  END DO
C
C
C Compute the TRIDIAGONAL operators, and the RHS source vector.
C
	  DO I=2,ND-1
	    TA(I)=-HL(I-1)
	    TC(I)=-HU(I)
	    TB(I)=DTAUONQ(I) + PSI(I) + HU(I-1) + HL(I)
	    VB(I)=-HS(I-1)
	    VC(I)=HS(I)
	    XM(I)=DTAUONQ(I)*SOURCE(I)*R(I)*R(I)
	  END DO
C
C Evaluate TA,TB,TC for boundary conditions
C
	  TC(1)=-F(2,ML)*Q(2)/DTAU(1)
	  TB(1)=F(1,ML)*Q(1)/DTAU(1) + HBC(1,ML) + PSI(1)
	  XM(1)=R(1)*R(1)*(HBC(2,ML)+SRCEBND)
	  TA(1)=0.0D0
	  VB(1)=0.0D0
	  VC(1)=0.0D0
C
	  TA(ND)=-F(ND-1,ML)*Q(ND-1)/DTAU(ND-1)
	  IF(DIF)THEN
	    TB(ND)=F(ND,ML)/DTAU(ND-1)
	    XM(ND)=DBB*R(ND)*R(ND)/3.0D0/TCHI(ND)
	  ELSE
	    TB(ND)=F(ND,ML)/DTAU(ND-1)+IN_HBC(ML)
	    XM(ND)=R(ND)*R(ND)*IC*(0.25D0+0.5D0*IN_HBC(ML))
	  END IF
	  TC(ND)=0.0D0
	  VB(ND)=0.0D0
	  VC(ND)=0.0D0
	  PSIPREV(ND)=0.0D0
C
	  IF(ML .NE. 1)THEN
	    XM(1)=XM(1) + PSIPREV(1)*JNUM1(1)
	    DO I=2,ND-1
	      XM(I)=XM(I) + VB(I)*HNUM1(I-1) + VC(I)*HNUM1(I)
	1          + PSIPREV(I)*JNUM1(I)
	    END DO
	    XM(ND)=XM(ND) + PSIPREV(ND)*JNUM1(ND)
	  END IF
C
C Solve for the radiation field along ray for this frequency.
C
	  CALL THOMAS(TA,TB,TC,XM,ND,1)
C
	  DO I=1,ND
	    JNU(I)=XM(I)
	  END DO
C
	  DO I=1,ND-1
	    HNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*HNUM1(I)
	  END DO
C 
C
C The J and H components of the radiation field have been found. We can thus
C begin the variation computation.
C
C Compute d{non-radiation field}/dchi matrix. NB - ML and NLF are
C required for the HBC and NBC vectors.
C
	  CALL EDDLINE_VAR(KI,RHS_dHdCHI,WORKMAT,
	1                SOURCE,TCHI,DTAU,R,SIGMA,
	1                F(1,ML),Q,HU,HL,HS,
	1                W,WPREV,PSI,PSIPREV,
	1                JNU,JNUM1,HNUM1,
	1                DEPTH_PRO,DBB,DIF,ML,NLF,ND,NM)
C
C Correct KI for variation of outer boundary condition with
C line source function.
C
	IF(ML .EQ. 1)THEN	!Prevent access of ML-1=0 term in HBC, and NBC.
	  T1=R(1)*R(1)*HBC(3,ML)
	ELSE
	  T1=R(1)*R(1)*(  HBC(3,ML)
	1                + GAM(1)*( HBC(3,ML)+SIGMA(1)*NBC(3,ML)
	1                    -HBC(3,ML-1)-SIGMA(1)*NBC(3,ML-1) )  )
	END IF
	KI(1,1,1)=KI(1,1,1)-T1*ETAL(1)/CHIL(1)/CHIL(1)
	KI(1,1,2)=KI(1,1,2)+T1/CHIL(1)
C
C Evaluate the intensity variations.
C                                  TX=d(u).d(chil,etal,chi,eta)
C                          and     TVX=d(v).d(chil,etal,chi,eta)
C
	  CALL UPTX_EDD(TX,TVX,KI,TA,TB,TC,
	1                  PSIPREV,VB,VC,ML_NE_ONE,ND,NM)
	  ML_NE_ONE=.TRUE.
C
C DEPTH_PHI must be equal to zero for the first frequency.
C
	  CALL UPTVX_EDD(TVX,TX,HU,HL,HS,
	1                   RHS_DHDCHI,DEPTH_PRO,ND,NM)
C
C Compute the "weights" to increment d{mean intensity}d{ , , , } arrays.
C Can't use DEPTH_PRO here as DEPTH_PRO is ZERO for ML=1. Must use
C PROF(1) for correct normalization.
C
	  DO I=1,ND
	    BAR_WTS(I)=LFQW(ML)*PROF(ML)
	  END DO
C
C Evaluate Diffusion variation. TXD is initially the value from the
C previous frequency.
C
	  IF(DIF)THEN
	    TXD(1)=PSIPREV(1)*TXD(1)
	    DO I=2,ND-1
	      TXD(I)= PSIPREV(I)*TXD(I)
	1             + VB(I)*TVXD(I-1) + VC(I)*TVXD(I)
	    END DO
	    TXD(ND)=R(ND)*R(ND)/3.0D0/TCHI(ND)+PSIPREV(ND)*TXD(ND)
C
C Solve for the radiation field along ray for this frequency.
C
	    CALL SIMPTH(TA,TB,TC,TXD,ND,1)
C
	    DO I=1,ND-1
	      TVXD(I)=HU(I)*TXD(I+1)-HL(I)*TXD(I)+HS(I)*TVXD(I)
	    END DO
C
	    DO I=1,ND
	      FQAFD(I)=FQAFD(I)+BAR_WTS(I)*TXD(I)
	    END DO
C
	  END IF	    	    !DIFF END
C
	  CALL MULT2D(FQAF,BAR_WTS,TX,ND,ND,NM)
C
C Save JNU and HNU for next frequency integration.
C
	  DO I=1,ND
	    JNUM1(I)=JNU(I)
	    HNUM1(I)=HNU(I)
	  END DO
C
	END DO
C
C 
C
C
C We have computed d(R^2 . J)/d{ }. Thus need to divide the variation matrix
C by r^2.
C
C
	DO I=1,ND
	  TA(I)=R(I)*R(I)
	  FQAFD(I)=FQAFD(I)/TA(I)
	END DO
	DO K=1,NM
	  DO J=1,ND
	    DO I=1,ND
	      FQAF(I,J,K)=FQAF(I,J,K)/TA(I)
	    END DO
	  END DO
	END DO
C
	IF(PROGDESC .NE. PROG_ID)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error - SCRATCH block corrupted in VAR_MOMJBAR'
	  STOP
	END IF
C
	RETURN
	END
