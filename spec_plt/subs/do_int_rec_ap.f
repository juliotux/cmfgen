	SUBROUTINE DO_INT_REC_AP(IP,P,YV,A,B,WIDTH,LENGTH,ZERO_YV,NP,NCF,NINS)
	IMPLICIT NONE
!
	INTEGER NP
	INTEGER NCF
	INTEGER NINS
	REAL*8 IP(NP,NCF)
	REAL*8 P(NP)
!
	REAL*8 A
	REAL*8 B
	REAL*8 WIDTH
	REAL*8 LENGTH
!
	REAL*4 YV(NCF)
	LOGICAL ZERO_YV
!
! Local variables
!
	REAL*8, ALLOCATABLE :: P_SHRT(:)
	REAL*8, ALLOCATABLE :: P_REV(:)
	REAL*8, ALLOCATABLE :: IP_REV(:)
	REAL*8, ALLOCATABLE :: dTHETA(:)
	REAL*8, ALLOCATABLE :: QW(:)
	REAL*8, ALLOCATABLE ::  INT_COEF(:)
	INTEGER, ALLOCATABLE :: INDX(:)

	INTEGER N_SP_P
	REAL*8 P_INS(10)
!
	REAL*8 P1,P2,P3,P4
!
	REAL*8 LEFT_EDGE
	REAL*8 RIGHT_EDGE
	REAL*8 BOT_EDGE
	REAL*8 TOP_EDGE
!
	REAL*8 PMIN
	REAL*8 PMAX
!
	REAL*8 T1,T2,T3
	REAL*8 TH1,TH2
	REAL*8 DELP
	REAL*8 PI
!
	INTEGER NP_SHRT
	INTEGER NP_REV
	INTEGER IMIN
	INTEGER IMAX
	INTEGER I,J,K,L,ML
	LOGICAL LINEAR
!
	LINEAR=.TRUE.
!	LINEAR=.FALSE.
!	NINS=(2*NINS+1)/2
	PI=ACOS(-1.0D0)
!
	IF(A .LT. 0.0D0 .OR. B .LT. 0.0D0)THEN
	  WRITE(6,*)'Error in DO_INT_REC_AP'
	  WRITE(6,*)'Aperture must be preferentially centered in top right quadrant'
	  STOP
	END IF
!
	IF(A .EQ. 0.0D0 .AND. B .EQ. 0.0D0)THEN
	ELSE IF(B-0.5D0*LENGTH .LT. 0.0D0)THEN
	  WRITE(6,*)'Error in DO_INT_REC_AP'
	  WRITE(6,*)'Aperture must be above X axis'
	  WRITE(6,*)'Call INT_RREC_AP first'
	  WRITE(6,*)'B=',B
	  WRITE(6,*)'LENGTH',0.5D0*LENGTH
	  STOP
	END IF
!	
	LEFT_EDGE=A-0.5D0*WIDTH
	RIGHT_EDGE=A+0.5D0*WIDTH
	BOT_EDGE=B-0.5D0*LENGTH
	TOP_EDGE=B+0.5D0*LENGTH
!
	WRITE(6,*)'Enterered DO_INT_REC_AP'
!
	IF(A .EQ. 0.0D0 .AND. B .EQ. 0.0D0)THEN
!
! Aperture is centered on the star.
!
	   PMIN=0.0D0
	   P1=0.5D0*MIN(WIDTH,LENGTH)
	   P2=0.5D0*MAX(WIDTH,LENGTH)
           PMAX=0.5D0*SQRT(WIDTH*WIDTH+LENGTH*LENGTH)
	   P_INS(1)=P1
	   P_INS(2)=P2
	   P_INS(3)=PMAX
	   N_SP_P=3
	   WRITE(6,*)'Listing of special P values'
	   WRITE(6,*)PMIN,PMAX,P1,P2,P(NP)
!
! Aperture is entirely within the right quadrant.
!
	ELSE IF(LEFT_EDGE .GE. 0.0D0 .AND. BOT_EDGE .GE. 0.0D0)THEN
	   P1=SQRT( LEFT_EDGE**2 + BOT_EDGE** 2)
	   P2=SQRT( LEFT_EDGE**2 + TOP_EDGE** 2)
	   P3=SQRT( RIGHT_EDGE**2 + TOP_EDGE** 2)
	   P4=SQRT( RIGHT_EDGE**2 + BOT_EDGE** 2)
!
	   PMIN=P1
	   PMAX=P3
	   P_INS(1)=PMIN
	   P_INS(2)=MIN(P2,P4)
	   P_INS(3)=MAX(P2,P4)
	   P_INS(4)=P3
	   N_SP_P=4
	   WRITE(6,*)'Listing of special P values'
	   WRITE(6,'(4ES16.8)')(P_INS(I),I=1,4)
!
! Aperture cross y-axis, but is in top half of plane.
! 
	ELSE IF(BOT_EDGE .GE. 0)THEN
	   P1=SQRT( LEFT_EDGE**2 + BOT_EDGE** 2)
	   P2=SQRT( LEFT_EDGE**2 + TOP_EDGE** 2)
	   P3=SQRT( RIGHT_EDGE**2 + TOP_EDGE** 2)
	   P4=SQRT( RIGHT_EDGE**2 + BOT_EDGE** 2)
!
	   PMIN=BOT_EDGE
	   PMAX=P3
	   P_INS(1)=PMIN
	   P_INS(2)=P1
	   P_INS(3)=TOP_EDGE
	   P_INS(4)=MIN(P2,P4)
	   P_INS(5)=MAX(P2,P4)
	   P_INS(6)=P3
	   IF(P_INS(3) .GT. P_INS(4))THEN
	     T1=P_INS(3)
	     P_INS(3)=P_INS(4)
             P_INS(4)=T1
	   END IF
	   N_SP_P=6
	   WRITE(6,*)'Listing of special P values'
	   WRITE(6,'(6ES16.8)')(P_INS(I),I=1,6)
!
	ELSE
	  WRITE(6,*)'Error in DO_INT_REC_AP'
	  WRITE(6,*)'Should never reach these instructions'
	  WRITE(6,*)'Chk programming'
	  STOP
	END IF
!
! Create a finer P grid to give higher accuracy.
!
	IMIN=1
	DO I=2,NP
	  IF(PMIN .LT. P(I))THEN
	    IMIN=I-1
	    EXIT
	  END IF
	END DO
!
	IMAX=NP-1
	DO I=2,NP
	  IF(PMAX .LT. P(I))THEN
	    IMAX=I-1
	    EXIT
	  END IF
	END DO
!
	IF(ALLOCATED(P_SHRT))DEALLOCATE(P_SHRT)
	ALLOCATE (P_SHRT(NP+10))
!
	K=1
	P_SHRT(1)=PMIN
	DO I=IMIN+1,IMAX
	  DO L=1,N_SP_P
	    T2=P(I)
	    IF(P_INS(L) .GT. P_SHRT(K) .AND. P_INS(L) .LT. T2)THEN
	      K=K+1
	      P_SHRT(K)=P_INS(L)
	    END IF
	  END DO
	  K=K+1
	  P_SHRT(K)=P(I)
	END DO
	IF(P_SHRT(K) .LT. PMAX)THEN
	  K=K+1
	  P_SHRT(K)=PMAX
	END IF
	NP_SHRT=K
!
	DO I=1,NP_SHRT
	  WRITE(71,*)I,P_SHRT(I)
	END DO
!
	IF(ALLOCATED(P_REV))DEALLOCATE(P_REV)
	ALLOCATE (P_REV(NP_SHRT*(NINS+1)+1))
	K=0
	DO I=1,NP_SHRT-1
	  DELP=(P_SHRT(I+1)-P_SHRT(I))/(NINS+1)
	  K=K+1
	  P_REV(K)=P_SHRT(I)
	  DO J=1,NINS
	    K=K+1
	    P_REV(K)=P_REV(K-1)+DELP
	  END DO
	END DO
	K=K+1
	P_REV(K)=P_SHRT(NP_SHRT)
	NP_REV=K
!
	WRITE(71,*)'P_REV'
	DO I=1,NP_REV
	  WRITE(71,*)I,P_REV(I)
	END DO
!
	ALLOCATE (QW(NP_REV))
	ALLOCATE (dTHETA(NP_REV))
!
	IF(A .EQ. 0.0D0 .AND. B .EQ. 0.0D0)THEN
	  I=1
	  DO WHILE(P_REV(I) .LE. P1)
	    dTHETA(I)=2.0D0*PI
	    IMIN=I
	    I=I+1
	  END DO
	  WRITE(6,*)'Done Theta 1'
!
	  I=IMIN+1
	  DO WHILE(P_REV(I) .LE. P2)
	    dTHETA(I)=4.0D0*(0.5D0*PI-ACOS(0.5D0*WIDTH/P_REV(I)))
	    IMIN=I
	    I=I+1
	  END DO
	  WRITE(6,*)'Done Theta 2'
!
	  I=IMIN+1
	  DO WHILE(I .LE. NP_REV)
	    dTHETA(I)=4.0D0*( ASIN(0.5D0*LENGTH/P_REV(I)) - ACOS(0.5D0*WIDTH/P_REV(I)) )
	    I=I+1
	  END DO
	  WRITE(6,*)'Done Theta 3'
	ELSE IF(LEFT_EDGE .GE. 0 .AND. BOT_EDGE .GE. 0)THEN
!
	  WRITE(6,*)'Beginning dTHETA comp.'
	  DO I=1,NP_REV
	    IF(P_REV(I) .EQ. 0)THEN
	      T1=0.0D0
	    ELSE IF(P_REV(I) .LE. P2)THEN
	      T1=ACOS(LEFT_EDGE/P_REV(I))
	    ELSE
	      T1=ASIN(TOP_EDGE/P_REV(I))
	    END IF
	    IF(P_REV(I) .EQ. 0)THEN
	      T2=0.0D0
	    ELSE IF(P_REV(I) .LE. P4)THEN
	      T2=ASIN(BOT_EDGE/P_REV(I))
	    ELSE
	      T2=ACOS(RIGHT_EDGE/P_REV(I))
	    END IF
	    dTHETA(I)=T1-T2
	  END DO
	  WRITE(6,*)'Done QW'
!
	ELSE IF(BOT_EDGE .GE. 0)THEN
!
	  WRITE(6,*)'Beginning dTHETA comp.'
	  DO I=1,NP_REV
	    IF(P_REV(I) .EQ. 0)THEN
	      T1=0.0D0
	    ELSE IF(P_REV(I) .LE. P1)THEN
	      T1=0.5D0*PI+ACOS(BOT_EDGE/P_REV(I))
	    ELSE IF(P_REV(I) .LE. P2)THEN
	      T1=0.5D0*PI+ASIN(ABS(LEFT_EDGE)/P_REV(I))
	    ELSE
	      T1=ASIN(ABS(TOP_EDGE)/P_REV(I))
	    END IF
	    IF(P_REV(I) .EQ. 0)THEN
	      T2=0.0D0
	    ELSE IF(P_REV(I) .LE. P4)THEN
	      T2=ASIN(BOT_EDGE/P_REV(I))
	    ELSE
	      T2=ACOS(RIGHT_EDGE/P_REV(I))
	    END IF
	    dTHETA(I)=T1-T2
	    IF(P_REV(I) .LT. P2 .AND. P_REV(I) .GT. TOP_EDGE)THEN
	       T1=ACOS(TOP_EDGE/P_REV(I))
	       dTHETA(I)=dTHETA(I)-2.0D0*T1
	    END IF
	  END DO
	END IF
!
	WRITE(71,*)' '
	WRITE(71,*)'dTHETA'
	DO I=1,NP_REV
	  WRITE(71,*)I,P_REV(I),dTHETA(I)
	END DO
!
	QW(:)=0.0D0
	IF(LINEAR)THEN
	  DO I=1,NP_REV-1
	    T1=0.5D0*(P_REV(I+1)-P_REV(I))
	    QW(I)=QW(I)+dTHETA(I)*P_REV(I)*T1
	    QW(I+1)=QW(I+1)+dTHETA(I+1)*P_REV(I+1)*T1
	  END DO
	ELSE
	  DO I=1,NP_REV-2,2
	    T1=(P_REV(I+1)-P_REV(I))/3.0D0
	    QW(I)=QW(I)+dTHETA(I)*P_REV(I)*T1
	    QW(I+1)=QW(I+1)+dTHETA(I+1)*P_REV(I+1)*T1*4.0D0
	    QW(I+2)=QW(I+2)+dTHETA(I+2)*P_REV(I+2)*T1
	  END DO
	END IF
!
	WRITE(71,*)' '
	WRITE(71,*)'QW'
	DO I=1,NP_REV
	  WRITE(71,'(I3,2X,3ES16.8)')I,P_REV(I),dTHETA(I),QW(I)
	END DO
!
	 WRITE(6,*)'Computed the quadrature weights'
!	
	IF(ZERO_YV)THEN
	  YV(1:NCF)=0.0D0	
	  WRITE(6,*)'Zeroed YV'
	END IF
	IF(ALLOCATED(IP_REV))DEALLOCATE(IP_REV)
	ALLOCATE (IP_REV(NP_REV))
!
	WRITE(6,*)'Beginning intepolation/integration section'
!
	IF(ALLOCATED(INDX))DEALLOCATE(INDX)
	ALLOCATE (INDX(NP_REV))
	IF(ALLOCATED(INT_COEF))DEALLOCATE(INT_COEF)
	ALLOCATE (INT_COEF(NP_REV))
!
	J=1
        DO I=1,NP_REV
500	  IF(P_REV(I) .LE. P(J+1))THEN
	    INDX(I)=J
	    INT_COEF(I)=(P_REV(I)-P(J))/(P(J+1)-P(J))
	  ELSE
	    J=J+1
	    GOTO 500
	  END IF
	END DO
!
	DO ML=1,NCF
!
! Perform simple linear interpolation.
!
	  J=2
          DO I=1,NP_REV
	    J=INDX(I)
	    T1=INT_COEF(I)
	    IP_REV(I)=(1.0D0-T1)*IP(J,ML)+T1*IP(J+1,ML)
	  END DO
!
!	WRITE(71,*)' '
!	WRITE(71,*)'IP_REV'
!	DO I=1,NP_REV
!	  WRITE(71,'(I3,2X,4ES16.8)')I,P_REV(I),dTHETA(I),QW(I),IP_REV(I)
!	END DO
!
!
! Do the quadrature.
!
	  IF(ML .EQ. 1)WRITE(71,'(A)')' Itegration sum'
	  DO I=1,NP_REV
	    YV(ML)=YV(ML)+QW(I)*IP_REV(I)
	    IF(ML .EQ. 1)WRITE(71,'(I3,2X,4ES16.8)')I,P_REV(I),dTHETA(I),QW(I),YV(ML)
	  END DO
	END DO
!
	RETURN
	END	
