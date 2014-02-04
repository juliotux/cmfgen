!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE VELOCITY_LAW(RVAL,ID,R,V,ND,VEL,BETA,DBETADR,GAMMA)
      IMPLICIT NONE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Subroutine to calculate the velocity, beta, dbetadr and gamma
! for the given radial point for a given velocity law.  The forms are:
!
!--------------------------------------------------------------------
!

      INTEGER ND
      INTEGER ID
!
      REAL*8 RVAL
      REAL*8 R(ND)
      REAL*8 V(ND)
!
      REAL*8 VEL
      REAL*8 BETA
      REAL*8 DBETADR
      REAL*8 GAMMA
!
      REAL*8, PARAMETER :: C_KMS=2.99792458D+05
      REAL*8 T1
!
!--------------------------------------------------------------------
!
! We use simple linear interpolaiton. Thus this procedure is valid for
! any velocity law --- not just a hubble expansion.
!
      T1=(RVAL-R(ID+1))/(R(ID)-R(ID+1))
      VEL=T1*V(ID)+(1.0D0-T1)*V(ID+1)
      BETA=VEL/C_KMS
      dBETAdR=(V(ID)-V(ID+1))/(R(ID)-R(ID+1))/C_KMS
      GAMMA=1.0D0/SQRT(1.0D0-BETA*BETA)
!
      RETURN
      END
