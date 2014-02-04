!
! Subroutine designed to replace the CMFGEN grid with the new grid.
! We use a subroutine to so we can access the grid in CMFGEN through MOD_CMFGEN.
!
	SUBROUTINE SET_NEW_GRID(REV_R,REV_V,REV_SIGMA,REV_ED,REV_CHI_ROSS,ND)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created 26-Jan-2007
!
	INTEGER ND
	REAL*8 REV_R(ND)
	REAL*8 REV_V(ND)
	REAL*8 REV_SIGMA(ND)
	REAL*8 REV_ED(ND)
	REAL*8 REV_CHI_ROSS(ND)
!
	R=REV_R
	V=REV_V
	SIGMA=REV_SIGMA
	ED=REV_ED
	ROSS_MEAN=REV_CHI_ROSS
!
	RETURN
	END
