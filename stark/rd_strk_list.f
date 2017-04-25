	MODULE MOD_STRK_LIST
	IMPLICIT NONE
!
! Altered 06-MAy-2014 : Put explicit stop in if using Lemke HI profile in IR
! Altered 24-Aug-2009 : Bug fix: LUER not being initialized for one error message.
! Altered 11-may-2004 : Bug fix: LST_V_PROF_LIMIT was not being sorted.
! Altered 03-Jan-2001 : LST_TYPE length increased from 10 to 12.
!
	INTEGER N_LST
!
	REAL*8, ALLOCATABLE :: LST_WAVE(:)
	REAL*8, ALLOCATABLE :: LST_GAM_RAD(:)
	REAL*8, ALLOCATABLE :: LST_GAM_COL(:)
	REAL*8, ALLOCATABLE :: LST_V_PROF_LIMIT(:)
!
	INTEGER, ALLOCATABLE :: LST_PID(:)
	INTEGER, ALLOCATABLE :: LST_NL(:)
	INTEGER, ALLOCATABLE :: LST_NUP(:)
!
	CHARACTER(LEN=12), ALLOCATABLE :: LST_TYPE(:)
	CHARACTER(LEN=5),  ALLOCATABLE :: LST_SPECIES(:)
!
	REAL*8, ALLOCATABLE :: INT_WRK(:)
	REAL*8, ALLOCATABLE :: DP_WRK(:)
	INTEGER, ALLOCATABLE :: VEC_INDX(:)
	CHARACTER*20, ALLOCATABLE :: CHAR_WRK(:)
!
	END MODULE MOD_STRK_LIST
!
	SUBROUTINE RD_STRK_LIST(LU)
	USE MOD_STRK_LIST
	IMPLICIT NONE
!
	INTEGER LU
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	INTEGER J,L
	INTEGER IOS
	INTEGER IPOS
	INTEGER LUER
	CHARACTER*80 STRING
!
	INTEGER ERROR_LU
	REAL*8 LAM_VAC
	EXTERNAL ERROR_LU,LAM_VAC
!
	OPEN(UNIT=LU,FILE='FULL_STRK_LIST',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error opening FULL_STRK_LIST in in RD_STRK_LIST'
	    WRITE(LUER,*)'This file specifies which stark tables to use.'
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  END IF
	  STRING=' '
	  DO WHILE( INDEX(STRING,'!Number of transitions') .EQ. 0)
	    READ(LU,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error reading number of transitions in RD_STRK_LIST'
	      WRITE(LUER,*)'IOSTAT=',IOS
	      STOP
	    END IF
	  END DO      
	  READ(STRING,*)N_LST
!
	  ALLOCATE (LST_WAVE(N_LST),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (LST_V_PROF_LIMIT(N_LST),STAT=IOS); LST_V_PROF_LIMIT=0.0D0
	  IF(IOS .EQ. 0)ALLOCATE (LST_GAM_RAD(N_LST),STAT=IOS); LST_GAM_RAD=0.0D0
	  IF(IOS .EQ. 0)ALLOCATE (LST_GAM_COL(N_LST),STAT=IOS); LST_GAM_RAD=0.0D0
	  IF(IOS .EQ. 0)ALLOCATE (LST_NL(N_LST),STAT=IOS);      LST_NL=0
	  IF(IOS .EQ. 0)ALLOCATE (LST_NUP(N_LST),STAT=IOS);     LST_NUP=0
	  IF(IOS .EQ. 0)ALLOCATE (LST_PID(N_LST),STAT=IOS);	LST_PID=0
	  IF(IOS .EQ. 0)ALLOCATE (LST_TYPE(N_LST),STAT=IOS);	LST_TYPE=' '
	  IF(IOS .EQ. 0)ALLOCATE (LST_SPECIES(N_LST),STAT=IOS);	LST_SPECIES=' '
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error allocating memory in RD_STRK_LIST'
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  END IF
!
	  L=1               
	  DO L=1,N_LST
10	    READ(LU,'(A)')STRING
	    IF(STRING(1:1) .EQ. '!')GOTO 10
	    STRING=ADJUSTL(STRING)
	    IPOS=INDEX(STRING,' ')
	    LST_SPECIES(L)=STRING(1:IPOS)
!
! Now get keywords.
!
	    STRING(1:)=STRING(IPOS+1:)
	    IPOS=INDEX(STRING,'[VW]=')
	    IF(IPOS .NE. 0)THEN
	      READ(STRING(IPOS+5:),*)LST_WAVE(L)
	    ELSE
	      IPOS=INDEX(STRING,'[AW]=')
	      IF(IPOS .EQ. 0)THEN
	        LUER=ERROR_LU()
	        WRITE(LUER,*)'Error in RD_STRK_LIST'
	        WRITE(LUER,*)'Wavelength data must be specified for all'
	        WRITE(LUER,*)'lines in FULL_STRK_LIST'
	        STOP
	      END IF
	      READ(STRING(IPOS+5:),*)LST_WAVE(L)
	      LST_WAVE(L)=LAM_VAC(LST_WAVE(L))
	    END IF
!
! Double brackets in ADJUSTL if for a bug in DEC F90 compiler.
!
	    IPOS=INDEX(STRING,'[TYPE]=')
	    IF(IPOS .EQ. 0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in RD_STRK_LIST'
	      WRITE(LUER,*)'Profile type must be specified for all'
	      WRITE(LUER,*)'lines in FULL_STRK_LIST'
	      STOP
	    END IF
	    IPOS=IPOS+7
	    STRING(IPOS:)=ADJUSTL( (STRING(IPOS:)) )
	    J=INDEX(STRING(IPOS:),' ')+IPOS-1
	    LST_TYPE(L)=STRING(IPOS:J-1)
!
	    IPOS=INDEX(STRING,'[PID]=')
	    IF(IPOS .NE. 0)READ(STRING(IPOS+6:),*)LST_PID(L)
!
! These values must be the same as in the STARK tables, but are
! not used to associate a PROFILE with a particular line.
!
	    IPOS=INDEX(STRING,'[NL]=')
	    IF(IPOS .NE. 0)READ(STRING(IPOS+5:),*)LST_NL(L)
	    IPOS=INDEX(STRING,'[NUP]=')
	    IF(IPOS .NE. 0)READ(STRING(IPOS+6:),*)LST_NUP(L)
!
	    IPOS=INDEX(STRING,'[G_COL]=')
	    IF(IPOS .NE. 0)READ(STRING(IPOS+8:),*)LST_GAM_COL(L)
!
	    IPOS=INDEX(STRING,'[G_RAD]=')
	    IF(IPOS .NE. 0)READ(STRING(IPOS+8:),*)LST_GAM_RAD(L)
!
! Use to set a fixed limit on the width of the intrinsic profile.
! If not present, the defaul vaule in SET_PROF_LIMITS is used.
! Value should be in km/s
!
	    IPOS=INDEX(STRING,'[V_LIM]=')
	    IF(IPOS .NE. 0)READ(STRING(IPOS+8:),*)LST_V_PROF_LIMIT(L)
!
	  END DO
!
	  IOS=0
	  DO WHILE(IOS .EQ. 0)   
	    READ(LU,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)EXIT
	    IF(STRING .NE. ' ' .AND. STRING(1:1) .NE. '!')THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Warning --- possible error in FULL_STRK_LIST'
	      WRITE(LUER,*)'Extra records in file'
	    END IF
	  END DO
	CLOSE(LU)
!
! Need to sort data into wavelength order for easier access.
!
	ALLOCATE (VEC_INDX(N_LST),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (INT_WRK(N_LST),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (DP_WRK(N_LST),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CHAR_WRK(N_LST),STAT=IOS)
	IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error allocating SORT arrays in RD_STRK_LIST'
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	END IF
!
	CALL INDEXX(N_LST,LST_WAVE,VEC_INDX,L_TRUE)
	CALL SORTDP(N_LST,LST_WAVE,VEC_INDX,DP_WRK)
	CALL SORTDP(N_LST, LST_V_PROF_LIMIT, VEC_INDX, DP_WRK)
	CALL SORTDP(N_LST,LST_GAM_RAD,VEC_INDX,DP_WRK)
	CALL SORTDP(N_LST,LST_GAM_COL,VEC_INDX,DP_WRK)
	CALL SORTINT(N_LST,LST_NL,VEC_INDX,INT_WRK)
	CALL SORTINT(N_LST,LST_NUP,VEC_INDX,INT_WRK)
	CALL SORTINT(N_LST,LST_PID,VEC_INDX,INT_WRK)
	CALL SORTCHAR(N_LST,LST_SPECIES,VEC_INDX,CHAR_WRK)
	CALL SORTCHAR(N_LST,LST_TYPE,VEC_INDX,CHAR_WRK)
!
	DO J=1,N_LST
	  IF(LST_WAVE(J) .GT. 1.3D+04 .AND. LST_TYPE(J) .EQ. 'LEMKE_HI' .AND. LST_NL(J) .EQ. 4)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)' '
	    WRITE(LUER,'(1X,80A)')('*',L=1,70)
	    WRITE(LUER,'(1X,80A)')('*',L=1,70)
	    WRITE(LUER,*)' '
	    WRITE(LUER,*)'Warning--- you should not use Lemke stark profiles for IR Bracket lines as they contain an error'
	    WRITE(LUER,*)'See: Repolust et al 2005, A&A 440, 261 (page 4)'
	    WRITE(LUER,*)' '
	    WRITE(LUER,'(1X,80A)')('*',L=1,70)
	    WRITE(LUER,'(1X,80A)')('*',L=1,70)
	    WRITE(LUER,*)' '
	    EXIT
	  END IF
	END DO
!
	DEALLOCATE (VEC_INDX)
	DEALLOCATE (INT_WRK)
	DEALLOCATE (DP_WRK)
	DEALLOCATE (CHAR_WRK)
!
	RETURN
	END
