C
C Subroutine to derive an interpolated temperature distribution from
C a previously found model distribution. The program allows for a change
C in luminosity and radius. The interpolation variable Y is defined by
C dY = chi*(R/r)**2 dr. Chi is the Rosseland mean optical depth, which is
C assumed to be due only to electron scattering. The electron density
C must be passed to the routine.
C
C Based on
C          J(Tau) = L / 16[piR]**2 { 3 int[0,Tau]  (R/r)^2 dTau + c }
C
C or
C            J(y) = L / 16[piR]**2 {y +c } where
C              dy = -chi*(R/r)**2 dr
C
C We adopt R as the radius at which Tau=1 (See Mihals, page 245).
C
      SUBROUTINE INIT_TEMP_V2(R,ED,CLUMP_FAC,T,LUM,
     &                          TAU_SWITCH,ND,LU,FILNAME)
	IMPLICIT NONE
C
C Altered 08-Jun-2002 - Bug fix in where we switch T correction.
C Altered 24-May-2001 - Bug fix. EDOLD(1:NDOLD) not (1:ND)
C Altered 19_mar-2001 - Bug fixed with Y interpolation scale. We now normalize
C                          the scale by the radius at which Tau=1.
C Altered 07-Jul-1997 - CLUMP_FAC_OLD is read in and used, if it is present.
C Altered 24-May-1996 - Work arrays made allocateable.
C Created 07-Apr-1989 - Based on TEMPDIST
C
	INTEGER ND,LU
	REAL*8 LUM               !Luminosity of star in Lsun
        REAL*8 R(ND)             !Radius Units 10^10 cm)
        REAL*8 ED(ND)            !Electron density
        REAL*8 CLUMP_FAC(ND)     !Volume filling factor
        REAL*8 T(ND)             !Temperature in units od 10^4 K
!
! Tau below which we (slowly) switch of the interpolation.
!
        REAL*8 TAU_SWITCH
	CHARACTER*(*) FILNAME
C
	REAL*8, ALLOCATABLE :: NEWED(:)
	REAL*8, ALLOCATABLE :: TAU(:)
C
        REAL*8, ALLOCATABLE :: ROLD(:)
	REAL*8, ALLOCATABLE :: EDOLD(:)
	REAL*8, ALLOCATABLE :: TOLD(:)
	REAL*8, ALLOCATABLE :: TAUOLD(:)
	REAL*8, ALLOCATABLE :: CLUMP_FAC_OLD(:)
!
! Local variables.
!
	INTEGER ERROR_LU,LU_ER
	EXTERNAL ERROR_LU
!
        INTEGER I,K,NOLD,NDOLD,IOS
        REAL*8 LUMOLD                  !Luminosity of old model in Lsun
        REAL*8 LDL,MLDL                !Used for interpolation.
        REAL*8 T1,DI,ION_FRAC,VEL      !Used when reading T_IN
        REAL*8 RTAU1		       !Radius in new model where TAU(e.s.)=1.0
        REAL*8 RTAU1_OLD               !Radius in old model where TAU(e.s.)=1.0
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
	  ALLOCATE (ROLD(I),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (EDOLD(I),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (TOLD(I),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (TAU(I),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (TAUOLD(I),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (NEWED(I),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CLUMP_FAC_OLD(I),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LU_ER=ERROR_LU()
	    WRITE(LU_ER,*)'Error in INIT_TEMP_V2'
	    WRITE(LU_ER,*)'Unable to allocate data arrays'
	    WRITE(LU_ER,*)'STAT=',IOS
	    STOP
	  END IF
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
C Determine radius at which optical depth is unity. We use NEWED as a
C temporary vector for ED*CLUMP_FAC.
C
        NEWED(1:ND)=ED(1:ND)*CLUMP_FAC(1:ND)
        TAU(1)=6.65D-15*NEWED(1)*R(1)
	K=0
	DO I=2,ND
          TAU(I)=TAU(I-1)+6.65D-15*(NEWED(I-1)+NEWED(I))*(R(I-1)-R(I))*0.5D0
          IF(TAU(I) .LE. 1.0D0)K=I
        END DO
	IF(K .EQ. 0)THEN
	  WRITE(6,*)'Error -- computed TAU is < 1 in INIT_TEMP_V2'
	  STOP
	END IF
        T1=(1.0D0-TAU(K))/(TAU(K+1)-TAU(K))
        RTAU1=T1*R(K+1)+(1.0D0-T1)*R(K)
C
C Compute the spherical optical depth scale. Assume atmosphere has constant
C V at outer boundary.
C
        DO I=1,ND
          NEWED(I)=NEWED(I)*RTAU1*RTAU1/R(I)/R(I)
	END DO
        TAU(1)=6.65D-15*NEWED(1)*R(1)/3.0D0
	DO I=2,ND
          TAU(I)=TAU(I-1)+6.65D-15*(NEWED(I-1)+NEWED(I))*(R(I-1)-R(I))*0.5D0
	END DO
C
C Determine radius at which optical depth is unity in old model. We first
C need to allow for clumping.
C
        EDOLD(1:NDOLD)=CLUMP_FAC_OLD(1:NDOLD)*EDOLD(1:NDOLD)
        TAUOLD(1)=6.65D-15*EDOLD(1)*ROLD(1)
	DO I=2,NDOLD
          TAUOLD(I)=TAUOLD(I-1)+6.65D-15*(EDOLD(I-1)+EDOLD(I))*
	1             (ROLD(I-1)-ROLD(I))*0.5D0
          IF(TAUOLD(I) .LE. 1.0D0)K=I
        END DO
        T1=(1.0D0-TAUOLD(K))/(TAUOLD(K+1)-TAUOLD(K))
        RTAU1_OLD=T1*ROLD(K+1)+(1.0D0-T1)*ROLD(K)
C
C Compute the spherical optical depth scale.
C
        DO I=1,NDOLD
          EDOLD(I)=EDOLD(I)*(RTAU1_OLD/ROLD(I))**2
        END DO
        TAUOLD(1)=6.65D-15*EDOLD(1)*ROLD(1)/3.0D0
	DO I=2,NDOLD
          TAUOLD(I)=TAUOLD(I-1)+6.65D-15*(EDOLD(I-1)+EDOLD(I))*
	1            (ROLD(I-1)-ROLD(I))*0.5D0
	END DO
C
C Interpolat the temperature onto the new grid, using the spherical optical
C depth scale.
C
C For Tau > TAU_SWITCH, we use the exact (approximate expression). For Tau < TAU_SWITCH, we
C reduce the correction, so that we keep the old temperature scale at small
C optical depths.
C
	K=2
        LDL=(LUM*RTAU1_OLD*RTAU1_OLD/LUMOLD/RTAU1/RTAU1)**0.25D0
	DO I=1,ND
	  IF(TAU(I) .LT. TAU_SWITCH)THEN
	    MLDL=1+(LDL-1)*TAU(I)/TAU_SWITCH
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
	    T(I)=MLDL*TOLD(NDOLD)*(TAU(I)/TAUOLD(NDOLD))**0.25D0
	  END IF
	END DO
C
C Free up storage.
C
	DEALLOCATE (ROLD)
	DEALLOCATE (EDOLD)
	DEALLOCATE (TOLD)
	DEALLOCATE (TAU)
	DEALLOCATE (TAUOLD)
	DEALLOCATE (NEWED)
	DEALLOCATE (CLUMP_FAC_OLD)
C
	CLOSE(UNIT=LU)
C
	RETURN
	END
