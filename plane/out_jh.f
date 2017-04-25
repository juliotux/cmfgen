	SUBROUTINE OUT_JH(RSQ_J,RSQ_H,H_INBC,H_OUTBC,NU,NCF,R,V,ND,INIT,OPTION)
	IMPLICIT NONE
!
! Altered: 04-May-2016.  R and V grid is now output at beginning when INIT is TRUE. 
!                         This ensures that R and V are current. Old values may have 
!                         been wrtten out with GREY option.  
	INTEGER ND
	INTEGER NCF
	REAL*8 NU
	REAL*8 R(ND)
	REAL*8 V(ND)
        REAL*8 RSQ_J(ND)
	REAL*8 RSQ_H(ND)
	REAL*8 H_INBC
	REAL*8 H_OUTBC
!
! Used to indicate that we are passing the first frequency.
!
	LOGICAL INIT
	CHARACTER(LEN=*), OPTIONAL :: OPTION
!
! Lcoal variables:
!
	REAL*8  T1
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
!
	INTEGER I		!Used as depth index
	INTEGER IOS		!I/O error identifier
	INTEGER REC_LENGTH 
	INTEGER LU_ER,ERROR_LU,WARNING_LU,LU_WARN
	EXTERNAL ERROR_LU,WARNING_LU
!
	REAL*8,  SAVE :: NU_STORE=0.0D0
	INTEGER, SAVE :: ST_IREC=6
	INTEGER, SAVE :: IREC=0
	INTEGER, SAVE :: LU_OUT=0
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
!
! The maximum number of elements in each record follows from:
!      ND   for RSQ_JNU
!      ND-1 for RSQ_HNU
!      3    for the two boundary conditions, and NU.
!
	LU_WARN=WARNING_LU()
	IF(FIRST_TIME)THEN
	  CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	  CALL GET_LU(LU_OUT,'GET_JH_AT_CURRENT_TIME_STEP')
          REC_LENGTH=WORD_SIZE*(2*ND+2)/UNIT_SIZE
          CALL WRITE_DIRECT_INFO_V3(ND,REC_LENGTH,'10-Jul-2006','JH_AT_CURRENT_TIME',LU_OUT)
          OPEN(UNIT=LU_OUT,FILE='JH_AT_CURRENT_TIME',STATUS='UNKNOWN',ACTION='WRITE',
	1        RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS)
          IF(IOS .NE. 0)THEN
	    LU_ER=ERROR_LU()
            WRITE(LU_ER,*)'Error opening JH_AT_CURRENT_TIME'
            WRITE(LU_ER,*)'IOS=',IOS
            STOP
	  END IF
	  DO IREC=1,ST_IREC-1
	    WRITE(LU_OUT,REC=IREC)0
	  END DO
	  WRITE(LU_OUT,REC=ST_IREC)(R(I),I=1,ND),(V(I),I=1,ND)
	  WRITE(LU_WARN,'(A)')' ST_IREC, IREC+2, NCF, ND, NU in OUT_JH'
	  WRITE(LU_WARN,*)ST_IREC,IREC,NCF,ND,NU
	  FIRST_TIME=.FALSE.
	END IF
!
	IF(PRESENT(OPTION))THEN
	  IF(OPTION .EQ. 'GREY')THEN
	    WRITE(LU_OUT,REC=ST_IREC+1)
	1          (RSQ_J(I),I=1,ND),(RSQ_H(I),I=1,ND-1),H_INBC,H_OUTBC
	    RETURN
	  ELSE IF(OPTION .EQ. 'NORMAL')THEN
	  ELSE
	    LU_ER=ERROR_LU()
	    WRITE(LU_ER,*)'Error in OUT_JH'
	    WRITE(LU_ER,*)'Invalid option passed'
	  END IF
	END IF
!
! Initialize indices. ST_IREC points to the record containing
! the first real data --- in this case R & V. IREC will contain
! the first record to be written for the frequency dependent J & H.
!
! The R,V write ensures R and V is uptodate if we have updated the R
! grid and used the GREY option.
!
	IF(INIT)THEN
	  WRITE(LU_OUT,REC=3)ST_IREC,NCF,ND
	  WRITE(LU_OUT,REC=ST_IREC)(R(I),I=1,ND),(V(I),I=1,ND)
	  IREC=ST_IREC+2   	!R,V, and frequency integrated J, H.
	  NU_STORE=0.0D0
	END IF
!
! Becasue we iterate, we may write the same frequency several times.
!
	IF(NU .EQ. NU_STORE)IREC=IREC-1
	WRITE(LU_OUT,REC=IREC)(RSQ_J(I),I=1,ND),(RSQ_H(I),I=1,ND-1),H_INBC,H_OUTBC,NU
	IREC=IREC+1
	NU_STORE=NU
!
	RETURN
	END
