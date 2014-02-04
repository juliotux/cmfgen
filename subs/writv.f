	SUBROUTINE WRITV(F,ND,A,LU)
	IMPLICIT NONE
!
! Altered 28-May-1996 : IMPLICIT NONE installed.
! Altered 30-APR-1985 : Now writes 10 columns instead of five across a page.)
!
	INTEGER ND,LU
	REAL*8 F(ND)
!
! Local variables.
!
	INTEGER I,J
	CHARACTER*(*) A
!
	WRITE(LU,'(//,1X,A,/)')A
	DO I=1,ND,10
	  WRITE(LU,'(1X,1P,10E12.4)')(F(J),J=I,MIN(I+9,ND),1)
	END DO
!
	RETURN
	END
