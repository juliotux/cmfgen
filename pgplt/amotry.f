      FUNCTION AMOTRY(P,Y,PSUM,MP,NP,NDIM,FUNK,IHI,FAC)
      IMPLICIT NONE
      INTEGER IHI,MP,NDIM,NP
      REAL*8 AMOTRY,FAC,P(MP,NP),PSUM(NP),Y(MP),FUNK
      INTEGER, PARAMETER :: NMAX=20
      EXTERNAL FUNK
!
! Uses function FUNK
!
      INTEGER J
      REAL*8 FAC1,FAC2,YTRY,PTRY(NMAX)
!
      FAC1=(1.D0-FAC)/NDIM
      FAC2=FAC1-FAC
      DO J=1,NDIM
        PTRY(J)=PSUM(J)*FAC1-P(IHI,J)*FAC2
      END DO
      YTRY=FUNK(PTRY)
!
      IF (YTRY .LT. Y(IHI)) THEN
        Y(IHI)=YTRY
        DO J=1,NDIM
          PSUM(J)=PSUM(J)-P(IHI,J)+PTRY(J)
          P(IHI,J)=PTRY(J)
	END DO
      END IF
      AMOTRY=YTRY
      RETURN
!
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 5
