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
	SUBROUTINE WIND_VEL_LAW_V2(R,V,SIGMA,VINF,BETA,BETA2,RMAX,R_TRANS,V_TRANS,
	1              dVdR_TRANS,VEL_TYPE,ND,ND_MAX)
	IMPLICIT NONE
!
! Altered 15-Jul-2007 -- Added VEL_TYPE=3, and now use LU_DIAG so lees ouput to OUTGEN.
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
	REAL*8 BETA2
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
	REAL*8 ALPHA
!
	REAL*8 T1,T2,T3
	INTEGER COUNT
	INTEGER I
	INTEGER LU_DIAG
	LOGICAL OUT_BOUNDARY
	LOGICAL FILE_OPENED
	LOGICAL VERBOSE
!
! Decide where diagnostic information will be written. HYDRO_ITERATION_INFO wil be open
! if this routine is called form DO_CMF_HYDRO_V2.
!
	LU_DIAG=6
	INQUIRE(FILE='HYDRO_ITERATION_INFO',NUMBER=I,OPENED=FILE_OPENED)
	IF(FILE_OPENED)LU_DIAG=I
	CALL GET_VERBOSE_INFO(VERBOSE)
	IF(VERBOSE)THEN
	  WRITE(LU_DIAG,*)' '
	  WRITE(LU_DIAG,*)'Called WIND_VEL_LAW_2 to set up the wind velocity'
	  WRITE(LU_DIAG,*)' '
	END IF
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
	  WRITE(LU_DIAG,*)'  Transition radius is',R_TRANS
	  WRITE(LU_DIAG,*)'Transition velocity is',V_TRANS
	  WRITE(LU_DIAG,*)'                 R0 is',RO
	  WRITE(LU_DIAG,*)'       Scale height is',SCALE_HEIGHT
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
!
	ELSE IF(VEL_TYPE .EQ. 2 .OR. VEL_TYPE .EQ. 3)THEN
	  IF(VERBOSE)THEN
	    WRITE(LU_DIAG,'(4X,A1,9(8X,A6))')'I','  R(I)','    T1','    T2','    T3',
	1                   'V(top)','V(bot)','dTOPdR','  dVdR','     V',' Sigma'
	  END IF
	  SCALE_HEIGHT = V_TRANS / (2.0D0 * DVDR_TRANS)
	  I=1
	  R(I)=R_TRANS
	  V(I)=V_TRANS
	  SIGMA(I)=R_TRANS*dVdR_TRANS/V_TRANS-1.0D0
	  ALPHA=2.0D0
	  IF(VEL_TYPE .EQ. 3)ALPHA=3.0D0
!
	  WRITE(LU_DIAG,*)'     Transition radius is',R_TRANS
	  WRITE(LU_DIAG,*)'   Transition velocity is',V_TRANS
	  WRITE(LU_DIAG,*)'                  dVdR is',dVdR_TRANS
	  WRITE(LU_DIAG,*)'                 SIGMA is',SIGMA(1)
	  WRITE(LU_DIAG,*)' Modified scale height is',SCALE_HEIGHT
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
	    T3=BETA+(BETA2-BETA)*T2
	    TOP = (VINF-ALPHA*V_TRANS) * T2**T3
	    BOT = 1.0D0 + (ALPHA-1.0D0)*exp( (R_TRANS-R(I))/SCALE_HEIGHT )
!
!NB: We drop a minus sign in dBOTdR, which is fixed in the next line.
!
	    dTOPdR = (VINF - ALPHA*V_TRANS) * BETA * T1 / R(I) * T2**(T3-1.0D0) +
	1                  T1*TOP*(BETA2-BETA)*(1.0D0+LOG(T2))/R(I)
	    dBOTdR=  (ALPHA-1.0D0)*exp( (R_TRANS-R(I))/SCALE_HEIGHT ) / SCALE_HEIGHT
!
	    TOP = ALPHA*V_TRANS + TOP
	    dVdR = dTOPdR / BOT  + TOP*dBOTdR/BOT/BOT
	    V(I) = TOP/BOT
            SIGMA(I)=R(I)*dVdR/V(I)-1.0D0
	    IF(VERBOSE)THEN
	      WRITE(LU_DIAG,'(I5,10ES14.4)')I,R(I),T1,T2,T3,TOP,BOT,dTOPdR,dVdR,V(I),SIGMA(I)
	    END IF
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
