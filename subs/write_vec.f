!
! This routine is essentially the same as WRITV, but it leaves less space.
! It cannot be called WRITV_V2, since a routine of that name already exists.
! I did not change WRITV to avoid possible problems in other peoples read
! routines.
!
	SUBROUTINE WRITE_VEC(F,ND,A,LU)
	IMPLICIT NONE
	INTEGER ND,LU
	REAL*8 F(ND)
!
! Local variables.
!
	INTEGER I,J
	CHARACTER*(*) A
!
	WRITE(LU,'(/,1X,A)')A
	DO I=1,ND,10
	  WRITE(LU,'(1X,1P,10E12.4)')(F(J),J=I,MIN(I+9,ND),1)
	END DO
!
	RETURN
	END
