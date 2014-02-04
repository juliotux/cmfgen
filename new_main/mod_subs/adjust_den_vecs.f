!
! Program to create a NEW_R_GRID which is equally spaced in LOG(Tau) where
! TAU is based on the FLUX mean opacity.
!
	SUBROUTINE ADJUST_DEN_VECS(R_OLD,ND)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created : 02-May-2004
!
	INTEGER ND
	REAL*8 R_OLD(ND)
	LOGICAL TRAPFORJ
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Local variables.
!
	REAL*8 LOG_R_OLD(ND)
	REAL*8 LOG_R(ND)
	REAL*8 TA(ND)
	REAL*8 TB(ND)
	INTEGER ISPEC
!
! We now need to regrid all the populations. All interpolations (except 
! sigma) are performed in the LOG-LOG plane. For SN this is ideal, since
! the density and velocity are power laws in r. 
!
	LOG_R_OLD=LOG(R_OLD)
	LOG_R=LOG(R)
!
! Now need to intepolate densities, and the cliumping factor, which all
! depend on the adopted R grid.
!
	TA(1:ND)=LOG(CLUMP_FAC(1:ND))
	CALL MON_INTERP(CLUMP_FAC,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	CLUMP_FAC(1:ND)=EXP(CLUMP_FAC(1:ND))
!
	TA(1:ND)=LOG(POP_ATOM(1:ND))
	CALL MON_INTERP(POP_ATOM,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	POP_ATOM(1:ND)=EXP(POP_ATOM(1:ND))
!
	TA(1:ND)=LOG(DENSITY(1:ND))
	CALL MON_INTERP(DENSITY,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	DENSITY(1:ND)=EXP(DENSITY(1:ND))
!
	DO ISPEC=1,NUM_SPECIES
	  IF(POP_SPECIES(1,ISPEC) .NE. 0.0D0)THEN
	    TA=LOG(POP_SPECIES(:,ISPEC))
	    CALL MON_INTERP(TB,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	    POP_SPECIES(:,ISPEC)=EXP(TB)
	  END IF
	END DO
!
	RETURN
	END
