C
C Routine to read in the logarithmic (base) 10 photoionization cross 
C sections for (n,l) and merged n states of hydrogen.
C
	MODULE HYD_BF_PHOT_DATA
	  INTEGER MAX_L_PQN
	  INTEGER N_PER_L
	  REAL*8 L_ST_U
	  REAL*8 L_DEL_U
	  REAL*8,    ALLOCATABLE ::  BF_L_CROSS(:)
	  INTEGER, ALLOCATABLE :: BF_L_INDX(:,:)
C
	  INTEGER MAX_N_PQN
	  INTEGER N_PER_N
	  REAL*8 N_ST_U
	  REAL*8 N_DEL_U
	  REAL*8,    ALLOCATABLE ::  BF_N_GAUNT(:)
	  INTEGER, ALLOCATABLE :: BF_N_INDX(:)
	END MODULE HYD_BF_PHOT_DATA
C
	SUBROUTINE RD_HYD_BF_DATA(LUIN,LUOUT,LUER)
	USE HYD_BF_PHOT_DATA
	IMPLICIT NONE
!
! Altered 22-Jun-2000 : Error ouput when HYD_L_DATA or GBF_N_DATA cannot be
!                        successfully opened.
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER LUIN,LUOUT,LUER
	INTEGER I,L,N,CNT,IOS
	INTEGER RD_L,RD_N
	CHARACTER*132 STRING
C
C Read in hydrogenic cross-section for (n,l) states.
C
        CALL GEN_ASCI_OPEN(LUIN,'HYD_L_DATA','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in RD_HYD_BF_DATA'
	  WRITE(LUER,*)'Unable to open HYD_L_DATA'
	  STOP
	END IF
C
C Read in header info.
C
	  STRING=' '
	  DO WHILE(INDEX(STRING,'Maximum principal quantum number') .EQ. 0)
	    READ(LUIN,'(A)')STRING
	    WRITE(LUOUT,'(A)')STRING
	  END DO
	  READ(STRING,*)MAX_L_PQN
C
	  READ(LUIN,'(A)')STRING
	  IF(INDEX(STRING,'Number of values') .EQ. 0)THEN
	    WRITE(LUER,*)'Error --- N_PER_L not found in RD_HYD_BF_DATA'
	    STOP
	  ELSE
	    READ(STRING,*)N_PER_L
	  END IF
C
	  READ(LUIN,'(A)')STRING
	  IF(INDEX(STRING,'L_ST_U') .EQ. 0)THEN
	    WRITE(LUER,*)'Error --- L_ST_U not found in RD_HYD_BF_DATA'
	    STOP
	  ELSE
	    READ(STRING,*)L_ST_U
	  END IF
C
	  READ(LUIN,'(A)')STRING
	  IF(INDEX(STRING,'L_DEL_U') .EQ. 0)THEN
	    WRITE(LUER,*)'Error --- L_DEL_U not found in RD_HYD_BF_DATA'
	    STOP                               
	  ELSE
	    READ(STRING,*)L_DEL_U
	  END IF
C
C Allocate necessary memorary.
C
	  ALLOCATE ( BF_L_CROSS(N_PER_L*MAX_L_PQN*(MAX_L_PQN+1)/2) )
	  ALLOCATE ( BF_L_INDX(MAX_L_PQN,0:MAX_L_PQN-1) )
C
	  CNT=0
	  DO N=1,MAX_L_PQN
	    DO L=0,N-1
	      READ(LUIN,*)RD_N,RD_L,I
	      IF(I .NE. N_PER_L)THEN
	        WRITE(LUER,*)'Error on RD_HYD_BF_DATA'
	        WRITE(LUER,*)'Invalid numer of elements for',N,L
	        STOP
	      END IF
	      IF(RD_N .NE. N .AND. RD_L .NE. L)THEN
	        WRITE(LUER,*)'Invalid N and L in RD_HYD_BF_DATA'
	        STOP
	      END IF
	      READ(LUIN,*)(BF_L_CROSS(CNT+I),I=1,N_PER_L)
	      BF_L_INDX(N,L)=CNT+1
	      CNT=CNT+N_PER_L
	    END DO
	  END DO
	CLOSE(LUIN)
C
C Now read in photoionization data for n levels.
C
        CALL GEN_ASCI_OPEN(LUIN,'GBF_N_DATA','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in RD_HYD_BF_DATA'
	  WRITE(LUER,*)'Unable to open GBF_N_DATA'
	  STOP
	END IF
C
C Read in header info.
C
	  STRING=' '
	  DO WHILE(INDEX(STRING,'Maximum principal quantum number') .EQ. 0)
	    READ(LUIN,'(A)')STRING
	    WRITE(LUOUT,'(A)')STRING
	  END DO
	  READ(STRING,*)MAX_N_PQN
C
	  READ(LUIN,'(A)')STRING
	  IF(INDEX(STRING,'Number of values') .EQ. 0)THEN
	    WRITE(LUER,*)'Error --- N_PER_N not found in RD_HYD_BF_DATA'
	    STOP
	  ELSE
	    READ(STRING,*)N_PER_N
	  END IF
C
	  READ(LUIN,'(A)')STRING
	  IF(INDEX(STRING,'N_ST_U') .EQ. 0)THEN
	    WRITE(LUER,*)'Error --- N_ST_U not found in RD_HYD_BF_DATA'
	    STOP
	  ELSE
	    READ(STRING,*)N_ST_U
	  END IF
C
	  READ(LUIN,'(A)')STRING
	  IF(INDEX(STRING,'N_DEL_U') .EQ. 0)THEN
	    WRITE(LUER,*)'Error --- N_DEL_U not found in RD_HYD_BF_DATA'
	    STOP                               
	  ELSE
	    READ(STRING,*)N_DEL_U
	  END IF
C
C Allocate necessary memorary.
C
	  ALLOCATE ( BF_N_GAUNT(N_PER_N*MAX_N_PQN) )
	  ALLOCATE ( BF_N_INDX(MAX_N_PQN) )
C
	  CNT=0
	  DO N=1,MAX_N_PQN
	    READ(LUIN,*)RD_N,I
	    IF(I .NE. N_PER_N)THEN
	      WRITE(LUER,*)'Error on RD_HYD_BF_DATA'
	      WRITE(LUER,*)'Invalid numer of elements for',N,L
	      STOP
	    END IF
	    IF(RD_N .NE. N)THEN
	      WRITE(LUER,*)'Invalid N and L in RD_HYD_PHOT_N'
	      STOP
	    END IF
	    READ(LUIN,*)(BF_N_GAUNT(CNT+I),I=1,N_PER_N)
	    BF_N_INDX(N)=CNT+1
	    CNT=CNT+N_PER_N
	  END DO
	CLOSE(LUIN)
C
	RETURN
	END
