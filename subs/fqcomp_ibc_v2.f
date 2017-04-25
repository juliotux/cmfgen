!
! This routine solves for the mean intensity as a function of depth using the
! Feautrier Technique. A Schuster or diffusion approaximation is used for the
! lower boundary condition. This routine may need to be in a loop so that
! the f values are iterated to convergence.
!
	SUBROUTINE FQCOMP_IBC_V2(TA,TB,TC,XM,DTAU,R,Z,P,NEWRJ,NEWRK,
	1            SOURCE,CHI,dCHIdr,AQW,AQW3,DBB,HBC_J,HBC_S,
	1            INBCNEW,IC,THK,INNER_BND_METH,NC,ND,NP,METHOD)
	IMPLICIT NONE
!
! Altered 15-Jan-2015: Added check to make sure there are no negative intensities.
! Altered 28-May-2013: J is now returned in XM.
! Altered 07-Jun-2010: Changed to THOMAS_RH to handle later stages of SN.
!                        Replaced DIF variable by INNER_BND_METH.
! Altered 26-May-2005 --- R*R-P*P replaved by (R-P)*(R+P)
! Altered 09-Dec-2001 --- PP replaced by P(LS)*P(LS)
! Altered 28-Oct-1997 --- Bug fix: COS corrected back to ACOS in TOR expression.
! Altered 24-May-1996 --- Call to DP_ZERO eliminated.
!                         IONE inserted in CALL to THOMAS
!                         Generic calls for EXP, COS
! Altered 12-Jun-1991 --- HBCNEW replaced by HBC_J and HBC_S. In order
!                         to improve convergence at outer boundary when wind
!                         is very thick. S1 deleted and use SOURCE(1) instead.
!                         Name changed from FQCOMP.
!
! Altered 21-Sep-1990 --- Bug fix: IBOUND was not being zeroed when THK
!                         option was off. Only effected direct computation of
!                         HBCNEW.
! Created 24-FEB-1987
!
	INTEGER NC,ND,NP,I,J,NI,LS
	REAL*8 TA(ND),TB(ND),TC(ND),XM(ND),R(ND),Z(ND)
	REAL*8 NEWRJ(ND),NEWRK(ND),DTAU(ND),dCHIdR(ND)
	REAL*8 SOURCE(ND),CHI(ND),AQW(ND,NP),AQW3(ND,NP),P(NP)
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER LU,LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	REAL*8 DBB,DBC,IBOUND,TOR,HBC_J,HBC_S,INBCNEW,IC,E1,E2,E3
	LOGICAL THK
	CHARACTER(LEN=6) METHOD
	CHARACTER(LEN=*) INNER_BND_METH
!
	 IF( INNER_BND_METH .NE. 'DIFFUSION' .AND.
	1   INNER_BND_METH .NE. 'SCHUSTER' .AND.
	1   INNER_BND_METH .NE. 'HOLLOW' .AND.
	1   INNER_BND_METH .NE. 'ZERO_FLUX' )THEN
	  WRITE(6,*)'Error in FQCOMP_IBC_V2'
	  WRITE(6,*)'Unrecognized INNER_BND_METH: ',TRIM(INNER_BND_METH)
	  STOP
	END IF
!
	CALL DERIVCHI(dCHIdr,CHI,R,ND,METHOD)
!
! Zero all parameters.
!
	NEWRJ(1:ND)=0.0D0
	NEWRK(1:ND)=0.0D0
	HBC_J=0.0D0
	HBC_S=0.0D0
	INBCNEW=0.0D0
!
! Enter loop for each impact parameter P.
!
	DO 2000 LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)THEN
	    NI=ND
	    DBC=DBB*DSQRT( (R(ND)-P(LS))*(R(ND)+P(LS)) )/R(ND)/CHI(ND)
	  END IF
!
	  IF(THK)THEN
	    IF(P(LS) .GT. 0.0D0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796D0-ACOS(P(LS)/R(1)))/P(LS)
	    ELSE
	      TOR=CHI(1)*R(1)
	    END IF
	    IBOUND=SOURCE(1)*(1.0D0-EXP(-TOR))
	  ELSE
	    IBOUND=0.0D0
	  END IF
!
! Compute Z and the optical depth scale DTAU for this imapct parameter.
!
	  IF(NI .GT. 1)THEN
	    DO I=1,NI
	      Z(I)=DSQRT( (R(I)-P(LS))*(R(I)+P(LS)) )
	    END DO
	    DO I=1,NI-1
	      DTAU(I)=0.5D0*(Z(I)-Z(I+1))*(CHI(I)+CHI(I+1)+(Z(I)-Z(I+1))
	1     *(dCHIdR(I+1)*Z(I+1)/R(I+1)-dCHIdR(I)*Z(I)/R(I))/6.0D0)
	    END DO
	  END IF
!
! Compute XM. Compute T ( a tridiagonal matrix) and store it as three vectors
! TA, TB and TC . This code is a combined version of XVECD and TCOMPD.
!
! To handle low otical depths, we use the Rybicki-Hummer formulation.
! In this formulation, TB(I)[tri matrix] = -[H(I)+TA(I)+TC(I)]
! We use TB as H in the following.
!
	  IF(NI .GT. 2)THEN
	    XM(1)=-IBOUND
	    TA(1)=0.0D0
	    TC(1)=1.0D0/DTAU(1)
	    TB(1)=1.0D0
	    DO I=2,NI-1
	      TA(I)=TC(I-1)
	      TC(I)=1.0D0/DTAU(I)
	      TB(I)=0.5D0*(DTAU(I-1)+DTAU(I))
	      XM(I)=-SOURCE(I)*(DTAU(I-1)+DTAU(I))*0.5D0
	    END DO
!
	    IF(LS .GT. NC)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-DTAU(NI-1)/2.0D0
	      XM(NI)=0.5D0*DTAU(NI-1)*SOURCE(NI)
	    ELSE IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=0.0D0
	      XM(NI)=DBC
	    ELSE IF(INNER_BND_METH .EQ. 'SCHUSTER')THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-1.0D0
	      XM(NI)=IC
	    ELSE
!
! HOLLOW core or ZERO_FLUX. For continnum, this is the same as for LS > NC. Done as separate option
! for clarity.
!
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-DTAU(NI-1)/2.0D0
	      XM(NI)=0.5D0*DTAU(NI-1)*SOURCE(NI)
	    END IF
	    TC(NI)=0.0D0
!
! Solve the tridiagonal system of equations.
!
	    CALL THOMAS_RH(TA,TB,TC,XM,NI,IONE)
!
	  ELSE IF(NI .EQ. 1)THEN
	    XM(1)=IBOUND
	  ELSE IF(NI .EQ. 2)THEN
	    E1=EXP(-DTAU(1))
	    E2=1.0D0-(1.0D0-E1)/DTAU(1)
	    E3=(1.0D0-E1)/DTAU(1)-E1
	    IF(DTAU(1) .LT. 1.0D-03)THEN
	      E2=DTAU(1)*0.5D0+DTAU(1)*DTAU(1)/6.0D0
	      E3=DTAU(1)*0.5D0-DTAU(1)*DTAU(1)/3.0D0
	    END IF
!
	    XM(2)=IBOUND*E1+SOURCE(2)*E2+SOURCE(1)*E3
	     XM(1)=0.5D0*(IBOUND+XM(2)*E1+SOURCE(1)*E2+SOURCE(2)*E3)
	  END IF
!
	  DO I=1,NI
	    IF(XM(I) .LT. 0.0D0)THEN
	      LUER=ERROR_LU()
	      CALL GET_LU(LU,'FQCOMP_IBC_V2')
	      WRITE(LUER,*)'Error in FQCOM_IBC_V2-- -ve intensities'
	      WRITE(LUER,*)'See unit LU'
	      WRITE(LU,*)'LS=',LS
	      WRITE(LU,'(6(10X,A))')'    XM','     Z','   CHI','dCHIdR','DTAU','SOURCE'
	      DTAU(NI)=0.0D0   !Set to zero as not used and not set.
	      DO J=1,NI
	        WRITE(LU,'(6ES16.8)')XM(J),Z(J),CHI(J),dCHIdR(J),DTAU(J),SOURCE(J)
	      END DO
	      STOP
	    END IF
	  END DO
!
! Update the FA and FB matrices (see notes).
!
	  DO I=1,NI
	    NEWRJ(I)=NEWRJ(I)+AQW(I,LS)*XM(I)
	    NEWRK(I)=NEWRK(I)+AQW3(I,LS)*XM(I)
	  END DO
!
	  HBC_J=HBC_J+AQW(1,LS)*XM(1)*Z(1)/R(1)
	  HBC_S=HBC_S+AQW(1,LS)*IBOUND*Z(1)/R(1)
	  IF(NI .EQ. ND)THEN
	    INBCNEW=INBCNEW + AQW(ND,LS)*
	1    (XM(ND)-(XM(ND)-XM(ND-1))*Z(ND)/R(ND)/DTAU(ND-1))
	  END IF
!
2000	CONTINUE
!
! Compute the factor for the outer boundary condition.
!
	HBC_J=HBC_J/NEWRJ(1)
	HBC_S=HBC_S/SOURCE(1)
	INBCNEW=INBCNEW/(2.0D0*NEWRJ(ND)-IC)
!
! Compute the new Feautrier factors.
!
	DO I=1,ND
	  NEWRK(I)=NEWRK(I)/NEWRJ(I)
	END DO
!
! Compute the Q factors from F. The Q values are stored in NEWRJ.
! We now compute Q in JFEAUNEW so that we nolonger need to save Q (only f).
! We still compute Q here because
!    (i) Q is passed (as NEWRJ vector) and hence is corrupted by this
!        routine. We still need Q, however, for the varaition
!        routines.
!   (ii) Q is now completely consistent with the updated f value.
!
	XM=NEWRJ
	CALL QFROMF(NEWRK,NEWRJ,R,TA,TB,ND)	!TA work vector
!
	RETURN
	END
