	SUBROUTINE OPEN_RW_EDDFACTOR(
	1     R,V,LANG_COORD,ND,
	1     R_EXT,V_EXT,LANG_COORD_EXT,ND_EXT,
	1     ACCESS_F,NEWMOD,COMPUTE_EDDFAC,USE_FIXED_J,FILENAME,LU_EDD)
	USE EDDFAC_REC_DEFS_MOD
	IMPLICIT NONE
!
! ACESS_F is the current record we are writing in EDDFACTOR.
! EDD_CONT_REC is the record in EDDFACTOR which points to the first
! record containing the continuum values.
!
	INTEGER ND,ND_EXT
	INTEGER ACCESS_F
	INTEGER LU_EDD
	LOGICAL NEWMOD
	LOGICAL COMPUTE_EDDFAC
	LOGICAL USE_FIXED_J
	CHARACTER(LEN=*) FILENAME
!
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 LANG_COORD(ND)
!
	REAL*8 R_EXT(ND)
	REAL*8 V_EXT(ND)
	REAL*8 LANG_COORD_EXT(ND)
!
	REAL*8 T1
!
! REC_SIZE     is the (maximum) record length in bytes.
! UNIT_SIZE    is the number of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the number of bytes used to represent the number.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
        INTEGER REC_SIZE
        INTEGER UNIT_SIZE
        INTEGER WORD_SIZE
        INTEGER N_PER_REC
        INTEGER RECORD_SIZE
!
	INTEGER IOS
	INTEGER I,K,J
	INTEGER IREC 
	INTEGER LUER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	INTEGER, PARAMETER :: IZERO=0
	CHARACTER(LEN=11) FILE_DATE
!
	ACCESS_F=0
	LUER=ERROR_LU()
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	RECORD_SIZE=WORD_SIZE*(ND_EXT+1)/UNIT_SIZE
!
! NB: If not ACCURATE, ND_EXT was set to ND. The +1 arises since we write
! NU on the same line as RJ. J is used to get the REC_LENGTH, while string
! will contain the date.
!
	IF(COMPUTE_EDDFAC)THEN
	ELSE
	  CALL READ_DIRECT_INFO_V3(K,J,FILE_DATE,FILENAME,LU_EDD,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error --- unable to open ',TRIM(FILENAME),'_INFO -- will compute new f'
	    COMPUTE_EDDFAC=.TRUE.
	    IOS=0
	  ELSE IF(.NOT. COMPUTE_EDDFAC .AND. K .NE. ND_EXT)THEN
	    WRITE(LUER,*)'Error with ',TRIM(FILENAME),'_INFO'
	    WRITE(LUER,*)'Incompatible number of depth points'
	    WRITE(LUER,*)'ND is',ND
	    WRITE(LUER,*)'ND_EXT is',ND_EXT
	    WRITE(LUER,*)'ND in ',TRIM(FILENAME),'_INFO is',K
	    WRITE(LUER,*)'You may need to delete the ',TRIM(FILENAME),' files'
	    STOP
	  END IF
	END IF
!
	IF(.NOT. COMPUTE_EDDFAC)THEN
	  OPEN(UNIT=LU_EDD,FILE=FILENAME,FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',RECL=RECORD_SIZE,IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    READ(LU_EDD,REC=FINISH_REC,IOSTAT=IOS)T1
	    IF(FILENAME .EQ. FILENAME .OR. FILENAME .EQ. 'ES_J_CONV')THEN
	      IF(T1 .EQ. 0.0D0 .OR. IOS .NE. 0)THEN
	        WRITE(LUER,'(/,A)')' Warning --- All EDDFACTOR values not'//
	1                      ' computed - will compute new F'
	        COMPUTE_EDDFAC=.TRUE.
	      END IF
	    END IF
	  ELSE
	    IF(.NOT. NEWMOD)THEN
	      WRITE(LUER,*)'Error opening ',TRIM(FILENAME),' - will compute new F'
	    END IF
	    COMPUTE_EDDFAC=.TRUE.
	  END IF
	END IF
!
	IF(COMPUTE_EDDFAC)THEN
	  IF(USE_FIXED_J)THEN
	    WRITE(LUER,*)'Error in CMFGEN_SUB'
	    WRITE(LUER,*)'Program will compute EDDFACTOR but this is'//
	1                      ' incompatable with US_FIXED_J=T'
	    STOP
	  END IF
!
          CALL WRITE_DIRECT_INFO_V3(ND_EXT,RECORD_SIZE,'23-Jan-2017',FILENAME,LU_EDD)
	  OPEN(UNIT=LU_EDD,FILE=FILENAME,FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='REPLACE',RECL=RECORD_SIZE)
	  WRITE(LU_EDD,REC=1)IZERO
	  WRITE(LU_EDD,REC=2)IZERO
	  WRITE(LU_EDD,REC=3)IZERO
	  WRITE(LU_EDD,REC=4)IZERO
!
! We set record 5 to zero, to signify that the Eddington factors are
! currently being computed. A non zero value signifies that all values
! have successfully been computed. (Consistent with old Eddfactor
! format a EDD_FAC can never be zero : Reason write a real number).
!
	  T1=0.0D0
	  WRITE(LU_EDD,REC=FINISH_REC)T1
	  ACCESS_F=INITIAL_ACCESS_REC
!
! ACCESS_REC will not be changed by CALL, since RV_REC has been updated.
!
	  IREC=INITIAL_RV_REC; WRITE(LU_EDD,REC=RV_REC)IREC
	  CALL OUT_RV_TO_EDDFACTOR(
	1        R,V,LANG_COORD,ND,
	1        R_EXT,V_EXT,LANG_COORD_EXT,ND_EXT,
	1        ACCESS_F,FILENAME,LU_EDD)
	END IF
!
	RETURN
	END
!
! We need to write out R, V and LANG_COORD at the end of each iteration since R may 
! have changed. This procedure will handle old format files, in which case R, V and
! the LANG_COORD will be writted at the end of the EDDFACTOR file.
!
	SUBROUTINE OUT_RV_TO_EDDFACTOR(
	1     R,V,LANG_COORD,ND,
	1     R_EXT,V_EXT,LANG_COORD_EXT,ND_EXT,
	1     ACCESS_F,FILENAME,LU_EDD)
	USE EDDFAC_REC_DEFS_MOD
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER ND_EXT
	INTEGER ACCESS_F
	INTEGER LU_EDD
!
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 LANG_COORD(ND)
!
	REAL*8 R_EXT(ND)
	REAL*8 V_EXT(ND)
	REAL*8 LANG_COORD_EXT(ND)
	CHARACTER(LEN=*) FILENAME
!
! REC_SIZE     is the (maximum) record length in bytes.
! UNIT_SIZE    is the number of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the number of bytes used to represent the number.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
        INTEGER REC_SIZE
        INTEGER UNIT_SIZE
        INTEGER WORD_SIZE
        INTEGER N_PER_REC
        INTEGER RECORD_SIZE
	INTEGER IREC
	INTEGER I
	INTEGER IOS
	INTEGER ERROR_LU
	INTEGER LUER
	LOGICAL TMP_LOG
!
	INQUIRE(UNIT=LU_EDD,OPENED=TMP_LOG)
	IF(.NOT. TMP_LOG)THEN
	  LUER=ERROR_LU()
	  CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	  I=WORD_SIZE*(ND_EXT+1)/UNIT_SIZE
	  OPEN(UNIT=LU_EDD,FILE=FILENAME,FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',RECL=I,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error oppeing in ',TRIM(FILENAME),' in OUT_RV_TO_EDDFACTOR'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	END IF
!
	READ(LU_EDD,REC=RV_REC)IREC
	IF(IREC .EQ. 0)THEN
	  IREC=ACCESS_F
	  WRITE(LU_EDD,REC=RV_REC)IREC
	  ACCESS_F=ACCESS_F+3
	END IF
!
	IF(ND_EXT .GT. ND)THEN
	  WRITE(LU_EDD,REC=IREC)R_EXT(1:ND_EXT)
	  WRITE(LU_EDD,REC=IREC+1)V_EXT(1:ND_EXT)
	  WRITE(LU_EDD,REC=IREC+2)LANG_COORD_EXT(1:ND_EXT)
	ELSE
	  WRITE(LU_EDD,REC=IREC)R(1:ND)
	  WRITE(LU_EDD,REC=IREC+1)V(1:ND)
	  WRITE(LU_EDD,REC=IREC+2)LANG_COORD(1:ND)
	END IF
!
	RETURN
	END
