!
! The simple routine computes:
!
!  Y=Y+X
!
! where Y and X are vectors of length NKT.
!
	SUBROUTINE PAR_VEC_SUM(Y,X,NKT)
	IMPLICIT NONE
	INTEGER NKT
	REAL*8 Y(NKT)
	REAL*8 X(NKT)
!
	INTEGER NOUT,NIN
	INTEGER J,IKT
!
	NOUT=16
        NIN=1+(NKT-1)/NOUT
!
!$OMP PARALLEL DO
        DO J=1,NOUT
          DO IKT=(J-1)*NIN+1,MIN(J*NIN,NKT)
            Y(IKT)=Y(IKT)+X(IKT)
          END DO
        END DO
!
	RETURN
	END
