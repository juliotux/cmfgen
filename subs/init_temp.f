C
C Subroutine to derive an interpolated temperature distribution from
C a previously found model distribution.  The program allows for a change
C in Luminosity and radius. Interpolation variable Y is defined by
C dY=chi*dr/r/r. Chi is the Rosseland mean optical depth, which is
C assumed to be due only to electron scattering. The electron density
C must be passed to the routine.
C
C The passed ED should be the average ED (i.e. ED*CLUMP_FAC)
C
	SUBROUTINE INIT_TEMP(R,ED,T,TAU,LUM,ND,LU,FILNAME)
	IMPLICIT NONE
C
C Altered 07-Jul-1997 - CLUMP_FAC is read in and used, if it is present.
C Altered 24-May-1996 - Work arrays made allocateable.
C Created 07-Apr-1989 - Based on TEMPDIST
C
	INTEGER ND,LU
	REAL*8 LUM
	REAL*8 R(ND),ED(ND),T(ND),TAU(ND)
	CHARACTER*(*) FILNAME
C
	REAL*8, ALLOCATABLE :: ROLD(:)
	REAL*8, ALLOCATABLE :: EDOLD(:)
	REAL*8, ALLOCATABLE :: TOLD(:)
	REAL*8, ALLOCATABLE :: TAUOLD(:)
	REAL*8, ALLOCATABLE :: NEWED(:)
	REAL*8, ALLOCATABLE :: CLUMP_FAC_OLD(:)
C
	INTEGER I,K,NOLD,NDOLD
	REAL*8 LUMOLD,LDL,MLDL
	REAL*8 T1,DI,ION_FRAC,VEL
	CHARACTER*132 STRING
	LOGICAL CLMP_PRES
C
	OPEN(UNIT=LU,FILE=FILNAME,STATUS='OLD')
C
C Check whether the clumping factor has also been written to file. The
C mere presence of a string contianung '!Format date' indicates that it has.
C
	  I=0
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	    I=I+1
	    READ(LU,'(A)')STRING
	  END DO
	  IF( INDEX(STRING,'!Format date') .NE. 0)THEN
	    CLMP_PRES=.TRUE.
	  ELSE
	    CLMP_PRES=.FALSE.
	    REWIND(LU)
	  END IF
C
	  READ(LU,*)T1,LUMOLD,NOLD,NDOLD
C
C Allocate required arrays.
C
	  I=MAX(ND,NDOLD)
	  ALLOCATE (ROLD(I))
	  ALLOCATE (EDOLD(I))
	  ALLOCATE (TOLD(I))
	  ALLOCATE (TAUOLD(I))
	  ALLOCATE (NEWED(I))
	  ALLOCATE (CLUMP_FAC_OLD(I))
C
	  IF(CLMP_PRES)THEN
	    DO I=1,NDOLD
	      READ(LU,*)ROLD(I),DI,EDOLD(I),TOLD(I),
	1                ION_FRAC,VEL,CLUMP_FAC_OLD(I)
	      READ(LU,*)(T1,K=1,NOLD)
	    END DO
	  ELSE
	    DO I=1,NDOLD
	      READ(LU,*)ROLD(I),DI,EDOLD(I),TOLD(I)
	      READ(LU,*)(T1,K=1,NOLD)
	      CLUMP_FAC_OLD(I)=1.0D0
	    END DO
	  END IF
	CLOSE(UNIT=LU)
C
	DO I=1,ND
	  NEWED(I)=ED(I)/R(I)/R(I)
	END DO
	TAU(1)=0.33333D0*NEWED(1)*R(1)*6.65D-15
	DO I=2,ND
	 TAU(I)=TAU(I-1)+(NEWED(I-1)+NEWED(I))*(R(I-1)-R(I))*0.5*6.65D-15
	END DO
C
	DO I=1,NDOLD
	  EDOLD(I)=CLUMP_FAC_OLD(I)*EDOLD(I)/ROLD(I)/ROLD(I)
	END DO
C
	TAUOLD(1)=0.33333D0*EDOLD(1)*ROLD(1)*6.65D-15
	DO I=2,NDOLD
	 TAUOLD(I)=TAUOLD(I-1)+(EDOLD(I-1)+EDOLD(I))
	1 *(ROLD(I-1)-ROLD(I))*0.5D0*6.65D-15
	END DO
C
	K=2
	LDL=(LUM/LUMOLD)**0.25
	T1=0.3333333D0/R(ND)/R(ND)
	DO I=1,ND
	  IF(TAU(I) .LT. T1)THEN
	    MLDL=1+(LDL-1)*TAU(I)/T1
	  ELSE
	    MLDL=LDL
	  END IF
	  IF(TAU(I) .LE. TAUOLD(1))THEN
	    T(I)=MLDL*TOLD(1)
	  ELSE IF(TAU(I) .LE. TAUOLD(NDOLD))THEN
10	    IF(TAU(I) .LE. TAUOLD(K))THEN
	      T(I)=MLDL*((TOLD(K)-TOLD(K-1))*(TAU(I)-TAUOLD(K-1))
	1     /(TAUOLD(K)-TAUOLD(K-1))+TOLD(K-1))
	    ELSE
	      K=K+1
	      GOTO 10
	    END IF
	  ELSE
	    T(I)=MLDL*TOLD(NDOLD)*(TAU(I)/TAUOLD(NDOLD))**0.25
	  END IF
	END DO
C
C Free up storage.
C
	DEALLOCATE (ROLD)
	DEALLOCATE (EDOLD)
	DEALLOCATE (TOLD)
	DEALLOCATE (TAUOLD)
	DEALLOCATE (NEWED)
	DEALLOCATE (CLUMP_FAC_OLD)
C
	CLOSE(UNIT=LU)
C
	RETURN
	END
