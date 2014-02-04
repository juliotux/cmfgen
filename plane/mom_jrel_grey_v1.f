!
! Routine to solve the frequency integrated, fully relativistic 
! transfer equations for a GREY atmosphere.
!
! Ref: Mihalas, ApJ, 237, 574
!      Solution of the Comoving-Frame Equation of Transfer in
!        Spherically Symmetric Flows. VI Relativistic flows.
!      See notes for implementation.
!
! The solution involves the use of 2 moment ratios. These moment
! ratios should be computed using a formal ray by ray solution. Routine
! should be called in a loop to converge the "MOMENT" ratios (or
! EDDINGTON) factors.
!
! Required MOMENT ratios:
!                	     F = K / J    (d=1,2,..., N)
!	                H_ON_J = H/J      (d=1,2,...,N)
!
! Also required are dBdR (=DBB) at the boundary. HBC should contain
! H/J at boundary.
!
! For the static solution, set V(1:ND)=0
!
! To neglect the Beta dJ/dR and Beta dHdR terms from the transfer
! equation (sometimes done for a supernovae model) set INCL_ADVEC_TERMSi
! to false.
!
	SUBROUTINE MOM_JREL_GREY_V1(ETA,CHI,ESEC,V,SIGMA,R,
	1		   H_ON_J,F,dlnJdlnR,
	1                  JNU,RSQHNU,HBC,IN_HBC,
	1                  DIF,DBB,IC,METHOD,COHERENT,
	1                  INCL_ADVEC_TERMS,INIT,ND)
	IMPLICIT NONE
!
! Created:   
!
	INTEGER ND
!
	REAL*8 R(ND)
	REAL*8 ETA(ND)
	REAL*8 CHI(ND)
	REAL*8 ESEC(ND)
	REAL*8 V(ND)
	REAL*8 SIGMA(ND)
!
! Moment ratio variables. All must be supplied.
!
	REAL*8 H_ON_J(ND)
	REAL*8 F(ND)
!
! If INIT is true, the initial vector values are irrelevant, and
! the variable is returned. If false, vector must be
! supplied, but it will be updated inside this routine.
!
	REAL*8 dlnJdlnR(ND)
!
! These values are computed, and returned.
!
	REAL*8 JNU(ND)
	REAL*8 RSQHNU(ND)
!
! Boundary conditions: Must be supplied.
!
	REAL*8 HBC			!H/J at outer boundary
	REAL*8 IN_HBC			!H/J at inner boundary
!
	REAL*8 DBB
	REAL*8 IC
	CHARACTER*6 METHOD
!
! COHERENT indicates whether the scattering is coherent. If it is, we
! explicitly take it into account. If COHERENT is FALSE, any electron
! scattering (or Rayleight scattering) term should be included directly 
! in the ETA that is passed to the routine.
!
	LOGICAL COHERENT
	LOGICAL DIF			!Diffusion approximation?
!
! If INCL_ADVEC_TERMS is true, the advection terms (dJ/dr and dH/dr)
!   are included. For supernova models, the ADVECTION terms are often 
!   neglected.
!
	LOGICAL INCL_ADVEC_TERMS
!
! If INIT is true, dlnJdlnR is initialized to zero. Otherwise the
!   passed values is used.
!
	LOGICAL INIT
!
! Local vectors
!
	REAL*8 TA(ND),TB(ND),TC(ND)
	REAL*8 CHI_H(ND),CHI_J(ND)
	REAL*8 DTAU_H(ND),DTAU_J(ND),DTAUONQ(ND)
	REAL*8 Q(ND),XM(ND),SOURCE(ND)
	REAL*8 HU(ND),HL(ND)
	REAL*8 BETA(ND),H_ADV_FAC(ND),GAM_REL(ND)
	REAL*8 GAM_REL_SQ(ND),CON_DELTA(ND)
	REAL*8 JOLD(ND),P_H(ND),P_J(ND)
	REAL*8 COH_VEC(ND)
!
! Local variables.
!
	REAL*8 T1,T2
	REAL*8 DAMP_FAC
	REAL*8 MAX_ER
	REAL*8 MAX_ER_LST_IT
	INTEGER LUER,ERROR_LU
	INTEGER I,J,IFAIL
	INTEGER IT_COUNT
	INTEGER, PARAMETER :: MAX_IT_COUNT=100
	INTEGER, SAVE :: LU_DIAG=0
!
	LOGICAL ACCURATE
	EXTERNAL ERROR_LU
!
! 
!
	IF(LU_DIAG .EQ. 0)THEN
	  CALL GET_LU(LU_DIAG,'MOM_JREL_GREY_V1')
	  OPEN(UNIT=LU_DIAG,FILE='GreyDiagnostics',STATUS='UNKNOWN')
	ELSE
	  CALL GET_LU(LU_DIAG,'MOM_JREL_GREY_V1')
	  OPEN(UNIT=LU_DIAG,FILE='GreyDiagnostics',STATUS='UNKNOWN',
	1              POSITION='APPEND')
	  WRITE(LU_DIAG,*)' '
	END IF
!
	LUER=ERROR_LU()
	DAMP_FAC=0.1D0
	MAX_ER_LST_IT=1.0D+10
!
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
	BETA(1:ND)=V(1:ND)/2.99794D+05
	H_ADV_FAC(1:ND)=BETA(1:ND)*H_ON_J(1:ND)
	GAM_REL_SQ(1:ND)=1.0D0/(1.0D0-BETA(1:ND)*BETA(1:ND))
	GAM_REL(1:ND)=SQRT(GAM_REL_SQ(1:ND))
	CON_DELTA(1:ND)=BETA(1:ND)/R(1:ND)
!
	IF(.NOT. INCL_ADVEC_TERMS)H_ADV_FAC(1:ND)=0.0D0
!
! 
!
! Zero relevant vectors and matrices.
!
	JNU(:)=0.0D0
	JOLD(:)=0.0D0
!
!*****************************************************************************
!
! CHI_H refers to the modified CHI term that multiplies H in the 1ST
! moment equation. We evaluate it at the grid points, and then use the
! averaging procedure used in the non-relativistic case.
!
! We divide CHI_H by an aditional GAM_REL since we substitute 
! gamma^2. r^2. H into the zero moment equation.
!
! CHI_J refers the the modified CHI term that multiplies J in the 0th
! moment equation.
!
	DO I=1,ND
	  CHI_H(I)=( CHI(I)/GAM_REL(I)+2.0D0*CON_DELTA(I)*(
	1       1.0D0+GAM_REL_SQ(I)*(SIGMA(I)+1) )  )/GAM_REL(I)
	  P_H(I)=1.0D0
	  CHI_J(I)=CHI(I)+GAM_REL(I)*CON_DELTA(I)*(
	1       3.0D0-F(I)+GAM_REL_SQ(I)*(1.0D0+SIGMA(I))*
	1                       (1.0D0+F(I)) )
	END DO
!
	SOURCE(1:ND)=ETA(1:ND)/CHI_J(1:ND)
	IF(COHERENT)THEN
	  COH_VEC(1:ND)=ESEC(1:ND)/CHI_J(1:ND)
	ELSE
	  COH_VEC(1:ND)=0.0D0
	END IF
!
! NB: We actually solve for gamma.r^2.J, not J.
!     We actually solve for (gamma.r)^2.H, not H.
!
! Compute the Q factors from F. The Q is used to make the J and K terms
! appearing in the first moment equation a perfect differential.
! Note that we integrate from the core outwards, and normalize Q to the
! value of unity at the core.
!
	DO I=1,ND
	  TA(ND-I+1)=( 3.0D0*F(I)-1.0D0-H_ADV_FAC(I)*(SIGMA(I)+1.0D0)+
	1      GAM_REL_SQ(I)*BETA(I)*BETA(I)*(SIGMA(I)+1.0D0)*
	1              (1.0D0-H_ADV_FAC(I)) )/(F(I)+H_ADV_FAC(I))/R(I)
	  TB(I)=R(ND-I+1)
	END DO
	CALL INTEGRATE(TB,TA,Q,IFAIL,ND)
!
	DO I=1,ND-1		!Q(ND) undefined after exiting INTEGRATE
	  TB(I)=EXP(Q(ND-I))
	END DO
!
! Scale by r^2
!
	DO I=1,ND-1
	  Q(I)=TB(I)/(R(I)/R(ND))**2
	END DO
	Q(ND)=1.0D0
!
! Compute optical depth scales.
!
	TA(1:ND)=CHI_H(1:ND)*Q(1:ND)
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU_H,TA,R,R,TB,ND)
!
	TA(1:ND)=CHI_J(1:ND)*Q(1:ND)
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU_J,TA,R,R,TB,ND)
!
	DO I=2,ND
	  DTAUONQ(I)=0.5D0*(DTAU_J(I)+DTAU_J(I-1))/Q(I)
	END DO
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=(F(I+1)+H_ADV_FAC(I+1))*Q(I+1)/P_H(I)/DTAU_H(I)
	  HL(I)=(F(I)+H_ADV_FAC(I))*Q(I)/P_H(I)/DTAU_H(I)
	END DO
!
! 
!
! As we don't know dlnJdlnR we need to iterate.
!
	ACCURATE=.FALSE.
!	IF(INIT)dlnJdlnR(1:ND)=-1.0D0
!
	DO IT_COUNT=1,MAX_IT_COUNT
!
	  IF(INCL_ADVEC_TERMS)THEN
	    P_J(1:ND)=1.0D0+CON_DELTA(1:ND)*
	1                  GAM_REL(1:ND)*dlnJdlnR(1:ND)/CHI_J(1:ND)
	  ELSE
	    P_J(1:ND)=1.0D0
	  END IF
!
! Compute the TRIDIAGONAL operators, and the RHS source vector. These
! vectors either depend  on dlnJdlnR etc, or are corrupted in the
! solution.
!
	  DO I=2,ND-1
	    TA(I)=-HL(I-1)
	    TC(I)=-HU(I)
	    TB(I)=DTAUONQ(I)*(P_J(I)-COH_VEC(I)) + HL(I) + HU(I-1)
	    XM(I)=DTAUONQ(I)*SOURCE(I)*R(I)*R(I)*GAM_REL(I)
	  END DO
!
! Evaluate TA,TB,TC for boundary conditions
!
	  TC(1)=-( F(2)+H_ADV_FAC(2))*Q(2)/DTAU_H(1)
	  TB(1)= ( F(1)+H_ADV_FAC(1))*Q(1)/DTAU_H(1) + HBC*P_H(1)*GAM_REL(1)
	  XM(1)=0.0D0
	  TA(1)=0.0D0
!
! NB: The expression for the flux is multiplied by GAM_REL_SQ since 
!     the first moment equation is written in terms of 
!     GAM_REL_SQ r^2 H which is required by the zero moment equation.
!
	  TA(ND)=-Q(ND-1)*(F(ND-1)+H_ADV_FAC(ND-1))/DTAU_H(ND-1)
	  IF(DIF)THEN
	    TB(ND)=(F(ND)+H_ADV_FAC(ND))/DTAU_H(ND-1)
	    XM(ND)=GAM_REL_SQ(ND)*DBB*R(ND)*R(ND)/3.0D0/CHI(ND)
	  ELSE
	    TB(ND)=(F(ND)+H_ADV_FAC(ND))/DTAU_H(ND-1)+IN_HBC*GAM_REL_SQ(ND)
	    XM(ND)=R(ND)*R(ND)*IC*(0.25D0+0.5D0*IN_HBC)*GAM_REL_SQ(ND)
	  END IF
	  TC(ND)=0.0D0
!
! Solve for the radiation field along ray for this frequency.
!
	  CALL THOMAS(TA,TB,TC,XM,ND,1)
!
! Check that no negative mean intensities have been computed.
!
	  J=0
	  DO I=1,ND
	    IF(XM(I) .LT. 0)THEN
	      J=J+1
	      XM(I)=ABS(XM(I))/10.0D0
	    END IF
	  END DO
	  IF(J .NE. 0)THEN
	    WRITE(LU_DIAG,*)'Error in MOM_JREL_GREY on iteration ',IT_COUNT
	    WRITE(LU_DIAG,*)'Zero intensity at ',J,'depths'
	   END IF
!
! Store J, correcting for the fact that we actually compute gamma r^2 J
!
	  DO I=1,ND
	    JNU(I)=XM(I)/R(I)/R(I)/GAM_REL(I)
	  END DO
!
! NB: RSQHNU contains r^2 H
!
	  DO I=1,ND-1
	    T1=0.5D0*(BETA(I)+BETA(I+1))
	    RSQHNU(I)=(HU(I)*XM(I+1)-HL(I)*XM(I))*(1.0D0-T1*T1)
	  END DO
	  RSQHNU(ND)=0.0D0
!
! Determine if J has converged. Not known if have included advections terms
! and are thus iterating on dlnJdlnR.
!
! To ensure convergence at very high velcities, we damp the changes.
! Since the estimate for dlnJdlnR may be very poor initially, we damp
! the changes very heavily for the first 3 iterations.
!
	  IF(.NOT. INCL_ADVEC_TERMS)THEN
	     ACCURATE=.TRUE.
	  ELSE
	    MAX_ER=0.0D0
	    DO I=1,ND
	      MAX_ER=MAX(MAX_ER,ABS((JOLD(I)-JNU(I))/JNU(I)))
	    END DO
	    WRITE(LU_DIAG,'(X,A,I3,A,ES14.6)')
	1                  'Maximum error in MOM_JREL_GREY_V1 on iteration ',
	1                  IT_COUNT,' is: ',MAX_ER
	    IF(MAX_ER .LT. 1.0D-06)ACCURATE=.TRUE.
	    JOLD(1:ND)=JNU(1:ND)
	    CALL DERIVCHI(TB,JNU,R,ND,'LINMON')
	    TB(1:ND)=R(1:ND)*TB(1:ND)/JNU(1:ND)
!	    IF(MAX_ER_LST_IT .GT. MAX_ER .AND. MAX_ER .LT. 0.001)THEN
!	      DAMP_FAC=MIN(0.8D0,DAMP_FAC+0.1D0)
!	    ELSE IF(MAX_ER_LST_IT .LT. MAX_ER)THEN          ! .AND. MAX_ER .GT. 0.05)THEN
!	      DAMP_FAC=MAX(0.1D0,DAMP_FAC-0.1D0)
!	    END IF
	    DAMP_FAC=0.2
	    dlnJdlnR(1:ND)=(1.0D0-DAMP_FAC)*dlnJdlnR(1:ND)+DAMP_FAC*TB(1:ND)
	    MAX_ER_LST_IT=MAX_ER
	  END IF
	  IF(ACCURATE)EXIT
	END DO
!
	IF(.NOT. ACCURATE)THEN
	  WRITE(LUER,*)'Error: too many iterations converging JGREY in MOM_JREL_GREY'
	  WRITE(LUER,*)'Number of iterations is',MAX_IT_COUNT
	  STOP
	END IF
!
! Check that "Radiative Equilibrium is being achived"'
!
	T1=HBC*GAM_REL(1)*XM(1)/R(1)/R(1)
	T2=GAM_REL_SQ(ND)*DBB/3.0D0/CHI(ND)
	CALL REGRID_H(TB,R,RSQHNU,T1,T2,ND,TC)
	TB=TB/GAM_REL					!g.r^2.H
	DO I=1,ND
!	  TA(I)=CON_DELTA(I)*( GAM_REL_SQ(I)*(SIGMA(I)+1.0D0)*
!	1                 (XM(I)+F(I)*XM(I)+BETA(I)*TB(I)) -
!	1                                  (SIGMA(I)+F(I))*XM(I) 
!	1                    )
	  TA(I)=CON_DELTA(I)*( XM(I)*(1.0D0-F(I)) +
	1           GAM_REL_SQ(I)*(SIGMA(I)+1.0D0)*(F(I)*XM(I)+BETA(I)*TB(I))
	1                    )
	  TC(I)=TB(I)+BETA(I)*XM(I)
	END DO
	CALL LUM_FROM_ETA(TA,R,ND)
	Q(1:ND)=0.0D0
	DO I=ND-1,1,-1
	  Q(I)=Q(I+1)+TA(I)
	END DO
!
	CLOSE(UNIT=LU_DIAG)
	OPEN(UNIT=LU_DIAG,STATUS='UNKNOWN',FILE='GREY_CHK')
	  WRITE(LU_DIAG,*)'HBC=',HBC
	  WRITE(LU_DIAG,'(4X,A,2X,8(3X,A),3X,A)')'I','        R','     Beta','  g.rsq.J',
	1                                            'gsq.rsq.H','  g.rsq.H','grsq.H+bJ',
	1                                            '    dInt','      Int','Cons Flux'
	  DO I=1,ND
	    WRITE(LU_DIAG,'(I5,2X,ES12.5,8ES12.3)')I,R(I),BETA(I),XM(I),RSQHNU(I),
	1                            TB(I),TC(I),TA(I),Q(I),TC(I)+Q(I)
	  END DO
	  T1=4.1274D-12*( TB(1)/GAM_REL(1) + GAM_REL(1)*BETA(1)*XM(1)*(1.0D0+F(1)) )
	  WRITE(LU_DIAG,'(A)')' '
	  WRITE(LU_DIAG,'(A,ES14.5)')'Observed flux computed by JREL_GREY_MOM_V1',T1
	  T1=4.1274D-12*TB(1)/GAM_REL(1)
	  WRITE(LU_DIAG,'(A,ES14.5)')'Comoving-frame flux computed by JREL_GREY_MOM_V1',T1
	  WRITE(LU_DIAG,'(A)')' '
	CLOSE(UNIT=LU_DIAG)
!
	RETURN
	END
