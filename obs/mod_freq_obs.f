!
! Placed large vectors (dimension with NCF_MAX and NLINE_MAX) in module MOD_FREQ_OBS.
!    These were being placed on the stack, and limiting the available stacks size.
!    Error in STACK found using -g debug option with -Mchkstk.
!
	MODULE MOD_FREQ_OBS
!
! Altered: 18-May-2015 : Added VEC_C4.
! Created: 22-Jan-2013
!
! NCF_MAX is the maximum number of continuum points (which will include line
! frequencies in blanketing mode) which can be treated. For small DOPPLER widths
! and large frequency ranges this might need to be increased.
!
	INTEGER, PARAMETER :: NCF_MAX=3000000
!
! Vectors for treating lines simultaneously with the continuum.
!
        INTEGER LINES_THIS_FREQ(NCF_MAX)
        INTEGER, ALLOCATABLE :: LINE_ST_INDX_IN_NU(:)
        INTEGER, ALLOCATABLE :: LINE_END_INDX_IN_NU(:)
!
        INTEGER, ALLOCATABLE :: LINE_LOC(:)     !Used to locate location of a particular line in the SIM vectors/arrays
!
! Continuum frequency variables and arrays.
!
	REAL*8 FQW(NCF_MAX)			!Frequency weights
        REAL*8 NU(NCF_MAX)              	!Continuum and line frequencies
        REAL*8 NU_EVAL_CONT(NCF_MAX)    	!Frequencies to evaluate continuum
        REAL*8 OBS(NCF_MAX)             	!Observers spectrum
!
        REAL*8 OBS_FREQ(NCF_MAX)                !Since N_OBS < NCF =< NCF_MAX
        REAL*8 OBS_FLUX(NCF_MAX)
!
! Arrays for performing LINE frequencies in numerical order
!
        INTEGER N_LINE_FREQ
        REAL*8, ALLOCATABLE :: VEC_FREQ(:)
        REAL*8, ALLOCATABLE :: VEC_STRT_FREQ(:)
        REAL*8, ALLOCATABLE :: VEC_OSCIL(:)
        REAL*8, ALLOCATABLE :: VEC_EINA(:)
        REAL*8, ALLOCATABLE :: VEC_ARAD(:)
        REAL*8, ALLOCATABLE :: VEC_C4(:)
        REAL*8, ALLOCATABLE :: VEC_DP_WRK(:)
        REAL*8, ALLOCATABLE :: VEC_VDOP_MIN(:)
!
        INTEGER, ALLOCATABLE :: VEC_INDX(:)
        INTEGER, ALLOCATABLE :: VEC_NL(:)
        INTEGER, ALLOCATABLE :: VEC_NUP(:)
        INTEGER, ALLOCATABLE :: VEC_MNL_F(:)
        INTEGER, ALLOCATABLE :: VEC_MNUP_F(:)
        INTEGER, ALLOCATABLE :: VEC_INT_WRK(:)
        INTEGER, ALLOCATABLE :: PROF_LIST_LOCATION(:)
!
        CHARACTER(LEN=6),  ALLOCATABLE :: VEC_SPEC(:)
        CHARACTER(LEN=6),  ALLOCATABLE :: VEC_TRANS_TYPE(:)
        CHARACTER(LEN=12), ALLOCATABLE :: PROF_TYPE(:)
        CHARACTER(LEN=12), ALLOCATABLE :: VEC_CHAR_WRK(:)
        CHARACTER(LEN=80), ALLOCATABLE :: VEC_TRANS_NAME(:)	!This routine is allocated in CMF_FLUX_SUB_V5 and then deallocated.
!
	END MODULE MOD_FREQ_OBS
!
	SUBROUTINE INIT_MOD_FREQ_OBS(NLINE_MAX)
	USE MOD_FREQ_OBS
	INTEGER NLINE_MAX
	INTEGER IOS
!
        ALLOCATE(VEC_FREQ(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_STRT_FREQ(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_OSCIL(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_EINA(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_ARAD(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_C4(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_DP_WRK(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_VDOP_MIN(NLINE_MAX),STAT=IOS)
!
        IF(IOS .EQ. 0)ALLOCATE(VEC_INDX(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_NL(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_NUP(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_MNL_F(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_MNUP_F(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_INT_WRK(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(PROF_LIST_LOCATION(NLINE_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(LINE_ST_INDX_IN_NU(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(LINE_END_INDX_IN_NU(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(LINE_LOC(NLINE_MAX),STAT=IOS)
!
        IF(IOS .EQ. 0)ALLOCATE(VEC_SPEC(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_TRANS_TYPE(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(PROF_TYPE(NLINE_MAX),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE(VEC_CHAR_WRK(NLINE_MAX),STAT=IOS)
!
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Error -- unable to allocate memory in INT_MOD_FREQ_OBS'
	  WRITE(6,*)'IOS=',IOS
	  STOP
	END IF
!
	RETURN
	END
