!
! This routine allocates the memory for the module used by SET_PROF.
! Data vectors are also initialized.
!
	SUBROUTINE INIT_PROF_MODULE(ND,NLINES,NFREQ)
	USE PROF_MOD
	IMPLICIT NONE
!
! Input: Note that NLINES only has to be large enough to handle overlapping
!        HI and HeII profiles.
!
	INTEGER ND		!Number of depth points
	INTEGER NLINES	!Maximum number of profiles to be stored
	INTEGER NFREQ		!Maximum number of frequencies / profile
!
	INTEGER IOS
	INTEGER ERROR_LU
	REAL*8 SPEED_OF_LIGHT
	REAL*8 FUN_PI
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
!
! Memory allocation. This should not be a problem, unless no more
! memory available. We simply insure all memory can be allocated.
!
	IOS=0
	ALLOCATE (PROF_STORE(ND,NFREQ,NLINES),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (NU_STORE(NFREQ,NLINES),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (AMASS_STORE(NLINES),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (NU_ZERO_STORE(NLINES),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (NL_STORE(NLINES),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (NUP_STORE(NLINES),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (NF_STORE(NLINES),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (LST_FREQ_LOC(NLINES),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (STORE_AVAIL(NLINES),STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE (PR_GRIEM(NFREQ),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (DWS_GRIEM(NFREQ),STAT=IOS)
!
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error allocating arrays in INIT_PROF_MOD'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
!
! Only vectors that needs to be initialized. Other locations will be
! set when line is available.
!
	NSTORE_MAX=NLINES
	NFREQ_MAX=NFREQ
	STORE_AVAIL(:)=.TRUE.
	PROF_STORE(:,:,:)=0.0D0
	LST_FREQ_LOC(:)=0
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	PI=FUN_PI()
	LUER=ERROR_LU()
!
	RETURN
	END
