!
! Subroutine to compute the wind velocity using an analytical formula,
! which is detrmined by VEL_TYPE. Only two options are currently implemented.
!
! A single parameter (in each velocity law) is adjusted so that the
! r, V, and dV/dR have specified values at the transition point.
!
! The routine returns R, V and SIGMA. The R grid is chose so that the spacing
! satisifies /\r < 0.3 r or /\v < 0.3 V.
!
	SUBROUTINE WIND_VEL_LAW_V1(R,V,SIGMA,VINF,BETA,RMAX,R_TRANS,V_TRANS,
	1              dVdR_TRANS,VEL_TYPE,ND,ND_MAX)
	IMPLICIT NONE
!
! Created 10-Aug-2006.
!
	INTEGER ND_MAX
	INTEGER ND
!
! Arrays output.
!
	REAL*8 R(ND_MAX)
	REAL*8 V(ND_MAX)
	REAL*8 SIGMA(ND_MAX)
!
! The following quantities must be specified on entry.
!
	REAL*8 VINF
	REAL*8 BETA
	REAL*8 RMAX
	REAL*8 R_TRANS
	REAL*8 V_TRANS
	REAL*8 dVdR_TRANS
	INTEGER VEL_TYPE
!
! Local variables.
!
	REAL*8 RO
	REAL*8 dR
	REAL*8 SCALE_HEIGHT
	REAL*8 TOP		!Numerator of velocity expression 
	REAL*8 BOT		!Denominator of velocity expression.
	REAL*8 dTOPdR,dBOTdR
	REAL*8 dVdR		!Velocoty gradient
!
	REAL*8 T1,T2
	INTEGER COUNT
	INTEGER I
	LOGICAL OUT_BOUNDARY
!
	OUT_BOUNDARY=.FALSE.
	COUNT=0
!
! We set the velocity grid from low velocities to high velocities. We put the
! grid into CMFGEN order at the end.
!
	IF(VEL_TYPE .EQ. 1)THEN
	  RO = R_TRANS * (1.0D0 - (2.0D0*V_TRANS/VINF)**(1.0D0/BETA) )
	  T1= R_TRANS * dVdR_TRANS / V_TRANS
	  SCALE_HEIGHT =  0.5D0*R_TRANS / (T1 - BETA*RO/(R_TRANS-RO) )
! 
	  WRITE(6,*)'  Transition radius is',R_TRANS
	  WRITE(6,*)'Transition velocity is',V_TRANS
	  WRITE(6,*)'                 R0 is',RO
	  WRITE(6,*)'       Scale height is',SCALE_HEIGHT
!
	  I=1
	  R(I)=R_TRANS
	  V(I)=V_TRANS
	  SIGMA(I)=R_TRANS*dVdR_TRANS/V_TRANS-1.0D0
	  DO WHILE (R(I) .LT. RMAX)
	    I=I+1
	    IF(I .GT. ND_MAX)THEN
	      WRITE(6,*)'Error in WIND_VEL_LAW: ND_MAX too small'
	      WRITE(6,*)'ND_MAX=',ND_MAX
	      STOP
	    END IF
	    IF(OUT_BOUNDARY)THEN
	      COUNT=COUNT+1
	      IF(COUNT .EQ. 1)THEN
	        R(I)=R(I-1)+0.6D0*dR
	      ELSE IF(COUNT .EQ. 2)THEN
	        R(I)=R(I-1)+0.27D0*dR
	      ELSE
	        R(I)=RMAX
	      END IF
	    ELSE 
	      T1=(SIGMA(I-1)+1.0D0)*V(I-1)/R(I-1)
	      dR=MIN(0.3D0*R(I-1),0.3D0*V(I-1)/T1)
	      R(I)=R(I-1)+dR
	      IF(R(I) + dR .GE. RMAX)THEN
	        OUT_BOUNDARY=.TRUE.
	        R(I)=(RMAX+R(I-1))/2.0D0
	        dR=RMAX-R(I)
	      END IF
	    END IF
!
            T1=RO/R(I)
            T2=1.0D0-T1
            TOP = VINF* (T2**BETA)
            BOT = 1.0D0 + exp( (R_TRANS-R(I))/SCALE_HEIGHT )
            V(I) = TOP/BOT
!                                                                                
!NB: We drop a minus sign in dBOTdR, which is fixed in the next line.
!                                                                               
            dTOPdR = VINF * BETA * T1 / R(I) * T2**(BETA - 1.0D0)
            dBOTdR=  exp( (R_TRANS-R(I))/SCALE_HEIGHT )  / SCALE_HEIGHT
            dVdR = dTOPdR / BOT  + V(I)*dBOTdR/BOT
            SIGMA(I)=R(I)*dVdR/V(I)-1.0D0
	  END DO
	ELSE IF(VEL_TYPE .EQ. 2)THEN
	  SCALE_HEIGHT = V_TRANS / (2.0D0 * DVDR_TRANS)
	  I=1
	  R(I)=R_TRANS
	  V(I)=V_TRANS
	  SIGMA(I)=R_TRANS*dVdR_TRANS/V_TRANS-1.0D0
!
	  WRITE(6,*)'  Transition radius is',R_TRANS
	  WRITE(6,*)'Transition velocity is',V_TRANS
	  WRITE(6,*)'               dVdR is',dVdR_TRANS
	  WRITE(6,*)'              SIGMA is',SIGMA(1)
	  WRITE(6,*)'       Scale height is',SCALE_HEIGHT
!
	  DO WHILE (R(I) .LT. RMAX)
	    I=I+1
	    IF(I .GT. ND_MAX)THEN
	      WRITE(6,*)'Error in WIND_VEL_LAW: ND_MAX too small'
	      WRITE(6,*)'ND_MAX=',ND_MAX
	      STOP
	    END IF
	    IF(OUT_BOUNDARY)THEN
	      COUNT=COUNT+1
	      IF(COUNT .EQ. 1)THEN
	        R(I)=R(I-1)+0.6D0*dR
	      ELSE IF(COUNT .EQ. 2)THEN
	        R(I)=R(I-1)+0.27D0*dR
	      ELSE
	        R(I)=RMAX
	      END IF
	    ELSE
	      T1=(SIGMA(I-1)+1.0D0)*V(I-1)/R(I-1)
	      dR=MIN(0.4D0*R(I-1),0.3D0*V(I-1)/T1)
	      R(I)=R(I-1)+dR
	      IF(R(I) + dR .GE. RMAX)THEN
	        OUT_BOUNDARY=.TRUE.
	        R(I)=(RMAX+R(I-1))/2.0D0
	        dR=RMAX-R(I)
	      END IF
	    END IF
!
	    T1=R_TRANS/R(I)
	    T2=1.0D0-T1
	    TOP = 2.0D0*V_TRANS + (VINF-2.0D0*V_TRANS) * T2**BETA
	    BOT = 1.0D0 + exp( (R_TRANS-R(I))/SCALE_HEIGHT )
	    V(I) = TOP/BOT
                                                                                
!NB: We drop a minus sign in dBOTdR, which is fixed in the next line.
                                                                                
	    dTOPdR = (VINF - 2.0D0*V_TRANS) * BETA * T1 / R(I) * T2**(BETA - 1.0D0)
	    dBOTdR=  exp( (R_TRANS-R(I))/SCALE_HEIGHT ) / SCALE_HEIGHT
	    dVdR = dTOPdR / BOT  + TOP*dBOTdR/BOT/BOT
            SIGMA(I)=R(I)*dVdR/V(I)-1.0D0
	  END DO
	ELSE
	  WRITE(6,*)'VEL_TYPE in WIND_VEL_LAW not recognized'
	  WRITE(6,*)'VEL_TYPE=',VEL_TYPE
	  STOP
	END IF
	ND=I
!
! Now reverse order of arrrays.
!
	DO I=1,ND/2
!
	  T1=R(I)
	  R(I)=R(ND-I+1)
	  R(ND-I+1)=T1
!
	  T1=V(I)
	  V(I)=V(ND-I+1)
	  V(ND-I+1)=T1
!
	  T1=SIGMA(I)
	  SIGMA(I)=SIGMA(ND-I+1)
	  SIGMA(ND-I+1)=T1
	END DO
!
	RETURN
	END
