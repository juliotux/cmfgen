	MODULE LINE_VEC_MOD
!
! Arrays for putting LINE frequencies in numerical order
!
	REAL*8,      ALLOCATABLE :: VEC_FREQ(:)             !NLINE_MAX
	REAL*8,      ALLOCATABLE :: VEC_OSCIL(:)            !NLINE_MAX
	REAL*8,      ALLOCATABLE :: VEC_EINA(:)             !NLINE_MAX
	REAL*8,      ALLOCATABLE :: VEC_DP_WRK(:)           !NLINE_MAX
	INTEGER,   ALLOCATABLE :: VEC_INDX(:)             !NLINE_MAX
	INTEGER,   ALLOCATABLE :: VEC_ID(:)               !NLINE_MAX
	INTEGER,   ALLOCATABLE :: VEC_NL(:)               !NLINE_MAX
	INTEGER,   ALLOCATABLE :: VEC_NUP(:)              !NLINE_MAX
	INTEGER,   ALLOCATABLE :: VEC_MNL_F(:)            !NLINE_MAX
	INTEGER,   ALLOCATABLE :: VEC_MNUP_F(:)           !NLINE_MAX
	INTEGER,   ALLOCATABLE :: VEC_INT_WRK(:)          !NLINE_MAX
	CHARACTER*6, ALLOCATABLE :: VEC_SPEC(:)             !NLINE_MAX
	CHARACTER*6, ALLOCATABLE :: VEC_CHAR_WRK(:)         !NLINE_MAX
	CHARACTER*6, ALLOCATABLE :: VEC_TRANS_TYPE(:)       !NLINE_MAX
!
	INTEGER, ALLOCATABLE :: LINE_ST_INDX_IN_NU(:)       !NLINE_MAX
	INTEGER, ALLOCATABLE :: LINE_END_INDX_IN_NU(:)      !NLINE_MAX
	INTEGER, ALLOCATABLE :: LINE_LOC(:)                !NLINE_MAX
!
! This will only be allocated temporarily.
!
	CHARACTER*80, ALLOCATABLE :: VEC_TRANS_NAME(:)
!
	END MODULE LINE_VEC_MOD