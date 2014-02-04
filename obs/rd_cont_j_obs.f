!
! Subroutine to return J interpolated from a file containing old J values.
! This routine is only to be used to compute J values during LAMBDA iterations
! when we have poor population estimates for some ionization stages/species.
! For use with CMF_FLUX only --- it is similar to RD_CONT_J used by CMFGEN.
!
! Altered 8-Mar-2005 : Bux fix.
!
	SUBROUTINE RD_CONT_J_OBS(RJ,FL,FREQ_INDX,FIRST_FREQ,
	1           LST_ITERATION,LUER,LU_EDD,ACCESS_F,ND)
	IMPLICIT NONE
!
! Altered 29-Jan-2010 : If frequency below minimum value, set J to J(last frequency).
! Created 28-Jan-2005
!
	INTEGER ND
	INTEGER ACCESS_F
	INTEGER FREQ_INDX
	INTEGER LU_EDD
	INTEGER LUER
!
	REAL*8 RJ(ND)
	REAL*8 FL
!
	LOGICAL LST_ITERATION
	LOGICAL FIRST_FREQ
!
! Local variables and arrays.
!
	REAL*8, SAVE :: LOW_FREQ, HIGH_FREQ
	REAL*8, ALLOCATABLE :: RJ_LOW(:)
	REAL*8, ALLOCATABLE :: RJ_HIGH(:)
	SAVE RJ_LOW
	SAVE RJ_HIGH
!
	REAL*8 T1
	INTEGER I
	INTEGER IOS
!
! Special treatment if first frequency. 
!
	IF(FIRST_FREQ)THEN
	  IF( .NOT. ALLOCATED(RJ_LOW) )THEN
	     ALLOCATE(RJ_LOW(1:ND))
	     ALLOCATE(RJ_HIGH(1:ND))
	  END IF
	  READ(LU_EDD,REC=ACCESS_F)(RJ_HIGH(I),I=1,ND),HIGH_FREQ
          ACCESS_F=ACCESS_F+1
	  READ(LU_EDD,REC=ACCESS_F)(RJ_LOW(I),I=1,ND),LOW_FREQ
          ACCESS_F=ACCESS_F+1
	END IF
!
! Get next frequency until we find the interpolation interval. Note
! that frequencies are ordered from highest to lowest.
!
	DO WHILE(FL .LT. LOW_FREQ)
	  HIGH_FREQ=LOW_FREQ
	  RJ_HIGH(1:ND)=RJ_LOW(1:ND)
	  READ(LU_EDD,REC=ACCESS_F,IOSTAT=IOS)(RJ_LOW(I),I=1,ND),LOW_FREQ
          IF(IOS. NE. 0)THEN
            LOW_FREQ=HIGH_FREQ
            RJ_LOW(1:ND)=RJ_HIGH(1:ND)
            WRITE(LUER,*)'Error reading EDDFACTOR in RD_CONT_J'
            WRITE(LUER,*)'LOW_FREQ=',LOW_FREQ
            WRITE(LUER,*)'FL=',FL
            WRITE(LUER,*)'ACCESS_F=',ACCESS_F
            EXIT
          END IF
          ACCESS_F=ACCESS_F+1
	END DO
!
! Do the interpolation: We use simple linear interpolation.
! If frequency higher than value in file, we adopt the first 
! file value.
!
        IF(FL .GE. HIGH_FREQ)THEN
	  RJ(1:ND)=RJ_HIGH(1:ND)
        ELSE IF(LOW_FREQ .EQ. HIGH_FREQ)THEN
	  RJ(1:ND)=RJ_LOW(1:ND)
	ELSE 
	  T1=(FL-LOW_FREQ)/(HIGH_FREQ-LOW_FREQ)
	  DO I=1,ND
	    RJ(I)=(1.0D0-T1)*RJ_LOW(I)+T1*RJ_HIGH(I)
	  END DO
	END IF
!
	RETURN
	END
