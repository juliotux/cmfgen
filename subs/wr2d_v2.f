C
C Routine to write a 2D Matrix out as a "MATRIX" . A maximum of ten
C numbers are written across the page.
C
	SUBROUTINE WR2D_V2(A,N,M,MES,SYMB,WR_INDEX,LU)
	IMPLICIT NONE
!
! Altered 10-Feb-2002: SYMB installed as pass option
!                      Index I & J now written out.
!                      Changed to V2
! Altered 13-Dec-1989 - Implicit none installed. I index written out.
!
	INTEGER N,M,LU
	REAL*8 A(N,M)
	LOGICAL WR_INDEX
	CHARACTER*(*) MES
	CHARACTER*(*) SYMB
C
	REAL*8 T1
	INTEGER MS,MF,ML,I,J
	CHARACTER*80 FORM
C
	WRITE(LU,'(/,1X,A)')MES
!
! If WR_INDEX is not set, we simply write the array values.
!
	IF(.NOT. WR_INDEX)THEN
	  MS=1
	  DO ML=0,M-1,10
	    MF=ML+10
	    IF(MF .GT. M)MF=M
	    WRITE(LU,'()')
	    DO I=1,N
	      WRITE(LU,'(1X,10ES12.4)')(A(I,J),J=MS,MF)
	    END DO
	    MS=MS+10
	  END DO
	  RETURN
	END IF
!
! Since WR_INDEX must be set, we need to determine the corect format to use.
!
! We Insert a 0 after the I so that FORM is always exactly 8 characters long
! before we enter the M section.
!
	IF(N .LT. 100)THEN
	   FORM='(1X,I02,'
	ELSE IF(N .LT. 1000)THEN
	   FORM='(1X,I03,'
	ELSE IF(N .LT. 10000)THEN
	   FORM='(1X,I04,'
	ELSE
	   T1=N
	   T1=LOG10(T1)+1
	   FORM='(1X,I'
	   WRITE(FORM(5:6),'(I2.2)')INT(T1)
	   FORM(8:8)=','
	END IF
!
	IF(M .LT. 10)THEN
	  FORM=TRIM(FORM)//'''('',I1,'')'','
	ELSE IF(M .LT. 100)THEN
	  FORM=TRIM(FORM)//'''('',I2,'')'','
	ELSE IF(M .LT. 1000)THEN
	  FORM=TRIM(FORM)//'''('',I3,'')'','
	ELSE IF(M .LT. 10000)THEN
	  FORM=TRIM(FORM)//'''('',I4,'')'','
	ELSE 
	   T1=M
	   T1=LOG10(T1)+1
	   FORM=TRIM(FORM)//'''('',I'
	   WRITE(FORM(14:15),'(I2.2)')INT(T1)
	   FORM=TRIM(FORM)//','')'','
	END IF
	IF(SYMB(1:1) .NE. ' ')THEN
	  FORM=TRIM(FORM)//''''//SYMB(1:1)//''',1X,10ES12.4)'
	ELSE
	  FORM=TRIM(FORM)//'1X,10ES12.4)'
	END IF
C
	MS=1
	DO ML=0,M-1,10
	  MF=ML+10
	  IF(MF .GT. M)MF=M
	  WRITE(LU,'()')
	  DO I=1,N
	    WRITE(LU,FORM)I,MS,(A(I,J),J=MS,MF)
	  END DO
	  MS=MS+10
	END DO
!
	RETURN
	END
