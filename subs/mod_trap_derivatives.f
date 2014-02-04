!
! Module to contain the data vectors with d(dCHIdR)/dCHI. Replaces common block
! TRAPDERIVATIVES.
!
! A(I)= d[dCHIdR(I)]/dCHI(I-1)
! B(I)= d[dCHIdR(I)]/dCHI(I)
! C(I)= d[dCHIdR(I)]/dCHI(I+1)
!
! Created: 02-Mar-1998. Allows arbitrarily large dimensions.
!
	MODULE MOD_TRAP_DERIVATIVES
	  INTEGER ND_TRAP
	  REAL*8, ALLOCATABLE :: A(:)
	  REAL*8, ALLOCATABLE :: B(:)
	  REAL*8, ALLOCATABLE :: C(:)
          DATA ND_TRAP/0/
	END MODULE MOD_TRAP_DERIVATIVES
