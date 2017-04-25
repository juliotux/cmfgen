	SUBROUTINE DIVIDE_BY_CONT(YV,NU,OBSF,NCF,NU_CONT,OBSF_CONT,NCF_CONT,LIN_INT)
	IMPLICIT NONE
!
	INTEGER NCF
	REAL*8 NU(NCF)
	REAL*8 OBSF(NCF)
	REAL*8 YV(NCF)
!
	INTEGER NCF_CONT
	REAL*8 NU_CONT(NCF_CONT)
	REAL*8 OBSF_CONT(NCF_CONT)
!
	LOGICAL LIN_INT
	LOGICAL UNEQUAL
	LOGICAL EQUAL
	EXTERNAL EQUAL
!
	REAL*8 T1,T2
	INTEGER I,J,L
!
	INTEGER, PARAMETER :: IONE=1
!	
	T1=1.0D-08
	I=1
	UNEQUAL=.FALSE.
	IF(NCF_CONT .NE. NCF)UNEQUAL=.TRUE.
	DO WHILE(.NOT. UNEQUAL .AND. I .LE. NCF_CONT)
	  IF( EQUAL(NU_CONT(I),NU(I),T1) )THEN
	    I=I+1
	  ELSE
	    UNEQUAL=.TRUE.
	  END IF
	END DO
!
	IF(UNEQUAL .AND. LIN_INT)THEN
	  L=1
	  DO I=1,NCF
	    IF(NU(I) .GT. NU_CONT(1))THEN
	      YV(I)=0.0
	    ELSE IF(NU(I) .LT. NU_CONT(NCF_CONT))THEN
	      YV(I)=0.0
	    ELSE 
	      DO WHILE (NU(I) .LT. NU_CONT(L+1))
	        L=L+1           
	      END DO
	      T1=(NU(I)-NU_CONT(L+1))/(NU_CONT(L)-NU_CONT(L+1))
	      T2=(1.0D0-T1)*OBSF_CONT(L+1)+T1*OBSF_CONT(L)
	      YV(I)=0.0
	      IF(T2 .NE. 0)THEN
	        T2=OBSF(I)/T2
	        IF(T2 .LT. 1.0E+020)YV(I)=T2
	      END IF
	    END IF
	  END DO
	ELSE IF(UNEQUAL)THEN
!
! We will use monotonic cubic interpolation. We first verify the range.
! I & J are temporary variables for the callt o MON_INTERP. I denotes the 
! first element. Initially J denotes the last element, then the numer of
! elements that can be interpolated.
!
	  I=1
	  DO WHILE(NU(I) .GT. NU_CONT(1))
	    I=I+1
	  END DO
	  J=NCF
	  DO WHILE(NU(J) .LE. NU_CONT(NCF_CONT))
	    J=J-1
	  END DO
	  J=J-I+1
!
	  YV(1:NCF)=0.0D0				!Temporary usage
	  CALL MON_INTERP(YV(I),J,IONE,NU(I),J,
	1            OBSF_CONT,NCF_CONT,NU_CONT,NCF_CONT)
	  DO I=1,NCF
	    IF(YV(I) .GT. 0)THEN
	      T2=OBSF(I)/YV(I)
	      IF(T2 .LT. 1.0E+30)YV(I)=T2
	    ELSE
	      YV(I)=0
	    END IF
	  END DO
	ELSE
	  DO I=1,NCF
	    YV(I)=0
	    IF(OBSF_CONT(I) .GT. 0)THEN
	      T2=OBSF(I)/OBSF_CONT(I)
	      IF(T2 .LT. 1.0E+20)YV(I)=T2
	    ELSE
	      YV(I)=0
	    END IF
	  END DO
	END IF
!
	RETURN
	END
