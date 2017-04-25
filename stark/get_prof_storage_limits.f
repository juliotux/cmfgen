	SUBROUTINE GET_PROFILE_STORAGE_LIMITS(NPROF_MAX,NFREQ_MAX,
	1       LINE_ST_INDX_IN_NU,LINE_END_INDX_IN_NU,PROF_TYPE,NLINES,NCF)
	IMPLICIT NONE
!
! These two variables are returned.
!
	INTEGER NPROF_MAX	!Maximum number of profiles to be stored at any one time
	INTEGER NFREQ_MAX	!Maximum number of frequencies for ANY profile.
!
! Required as input:
!
	INTEGER NLINES
	INTEGER NCF
	INTEGER LINE_ST_INDX_IN_NU(NLINES)
	INTEGER LINE_END_INDX_IN_NU(NLINES)
	CHARACTER(LEN=*) PROF_TYPE(NLINES)
!
! Local variables
!
	INTEGER K
	INTEGER NL, NL_SAVE
	INTEGER ML
	INTEGER LN_INDX
	INTEGER NX
!
! Get maximum number of frequencies required by any intrinsic line profile.
!
	NFREQ_MAX=0.0D0
	DO NL=1,NLINES
	  IF(PROF_TYPE(NL)(1:3) .NE. 'DOP' .AND. PROF_TYPE(NL) .NE. 'VOIGT')
	1      NFREQ_MAX=MAX(NFREQ_MAX,LINE_END_INDX_IN_NU(NL)-LINE_ST_INDX_IN_NU(NL)+1)
	END DO
!
! Get the maximum number of lines that overlap, and that require storage of their intrinsic profile.
! Line computed with Doppler or Voigt profiles do not require storage.
!
	NPROF_MAX=0
	LN_INDX=1
	DO ML=1,NCF
	  NX=0
	  K=LN_INDX
	  DO NL=K,NLINES
	    NL_SAVE=NL
	    IF(PROF_TYPE(NL)(1:3) .NE. 'DOP' .AND. PROF_TYPE(NL) .NE. 'VOIGT' .AND.
	1      ML .GE. LINE_ST_INDX_IN_NU(NL) .AND.
	1      ML .LE. LINE_END_INDX_IN_NU(NL) )NX=NX+1
	    IF(LINE_END_INDX_IN_NU(NL) .LE. ML .AND. NL .EQ. LN_INDX)LN_INDX=LN_INDX+1
	    IF(LINE_ST_INDX_IN_NU(NL) .GT. ML)EXIT
	  END DO
	  NPROF_MAX=MAX(NX,NPROF_MAX)
	  IF(LN_INDX .EQ. NLINES)EXIT
	END DO
!
	IF(NPROF_MAX .NE. 0 .OR. NFREQ_MAX .NE. 0)THEN
	  WRITE(6,*)'Maximum number of profiles to be stored is ',NPROF_MAX
	  WRITE(6,*)'Maximum number of frequncies per profile is',NFREQ_MAX
	END IF
!
	RETURN
	END
