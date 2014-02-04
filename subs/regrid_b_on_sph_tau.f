!
! Program reads in departire coefficient estimates from a file and uses
! linear interpolation in the log-tau plane to lay these estimates
! on the grid. The electron density is assumed to be known on entering 
! the program. The ion density is also returned.
!
! The subroutine uses the same tau scale as in INIT_TEMP_V2 for the
! interpolation. Thus the T's and d.c's should be obtained from the same
! locations.
!
! Routine also has the advantage over REGRID_B_ON_NE in that it treats
! non-monotonic electron distributions.
!
	SUBROUTINE REGRID_B_ON_SPH_TAU(DHEN,ION,ED,CLUMP_FAC,R,N,ND,LU,FILNAME)
	IMPLICIT NONE
!
! Created 07-Jun-2005 - Based on REGRID_B_ON_NE
!
	INTEGER N,ND,LU
	REAL*8 DHEN(N,ND)
	REAL*8 ION(ND)
	REAL*8 CLUMP_FAC(ND)
	REAL*8 ED(ND)
	REAL*8 R(ND)
	CHARACTER(LEN=*) FILNAME
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	REAL*8 NEW_ED(ND)
	REAL*8 TAU(ND)
!
! Old values read in from file. HYD is the old d.c.'s
!
	REAL*8, ALLOCATABLE :: HYD(:,:)
	REAL*8, ALLOCATABLE :: ROLD(:)
	REAL*8, ALLOCATABLE :: OLDED(:)
	REAL*8, ALLOCATABLE :: OLDTAU(:)
	REAL*8, ALLOCATABLE :: OLDION(:)
!
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
!
	INTEGER I,J,K
	INTEGER NX,NXST,NX_END,NZ
	INTEGER NOLD,NDOLD,IOS
	REAL*8 CONST,TSTAROLD,RPOLD
	REAL*8 ION_FRAC,VEL,OLD_CLUMP_FAC
	REAL*8 RTAU1,RTAU1_OLD,T1
	LOGICAL CLMP_PRES
	LOGICAL CHECK_DC
	CHARACTER*80 STRING
!
! Read in values from previous model.
!
	OPEN(UNIT=LU,STATUS='OLD',FILE=FILNAME,ACTION='READ',IOSTAT=IOS)
        IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error opening ',TRIM(FILNAME),' in REGRID_B_ON_NE'
          WRITE(LUER,*)'IOS=',IOS
          STOP
        END IF
!
! Check whether the clumping factor has also been written to file. The
! mere presence of a string containing '!Format date' indicates that it has.
!
	I=0
	STRING=' '
	DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	  I=I+1
	  READ(LU,'(A)')STRING
	END DO
	IF( INDEX(STRING,'!Format date') .NE. 0)THEN
	  CLMP_PRES=.TRUE.
	  CHECK_DC=.FALSE.
	ELSE
	  CHECK_DC=.TRUE.
	  CLMP_PRES=.FALSE.
	  REWIND(LU)
	END IF

	READ(LU,*)RPOLD,TSTAROLD,NOLD,NDOLD
!
! Allocate necessary memory.
!
	I=MAX(NOLD,NDOLD,N,ND)
	ALLOCATE (HYD(NOLD,NDOLD))
	ALLOCATE (ROLD(I))
	ALLOCATE (OLDED(I))
	ALLOCATE (OLDTAU(I))
	ALLOCATE (OLDION(I))
	ALLOCATE (TA(I))
	ALLOCATE (TB(I))
!
! NZ defines the number of atomic levels which can be found by
! direct interpolation
!
	NZ=N
	IF(NZ .GT. NOLD)NZ=NOLD
!
! TA is used for both R and T but note that they
! are not used and hence may be overwritten.
!
	DO I=1,NDOLD
	  IF(CLMP_PRES)THEN
	    READ(LU,*)ROLD(I),OLDION(I),OLDED(I),TA(I),ION_FRAC,VEL,OLD_CLUMP_FAC
	    OLDED(I)=OLDED(I)*OLD_CLUMP_FAC
	  ELSE
	    READ(LU,*)ROLD(I),OLDION(I),OLDED(I),TA(I)
	    OLDED(I)=OLDED(I)
	  END IF
	  READ(LU,*)(TA(J),J=1,NOLD)
	  DO J=1,NOLD
	    HYD(J,I)=TA(J)
	  END DO
	  OLDION(I)=LOG(OLDION(I))
	END DO
!
! Check to see if departure coefficients, or one minus departure
! coefficients. Assumes that at depth, the highest level should
! have b=1.
!
	IF(.NOT. CHECK_DC)THEN
	  CONST=0.0D0
	ELSE IF( ABS((HYD(NZ,NDOLD)-1.0D0)) .LT. 0.1)THEN
	  CONST=0.0D0
	ELSE IF( ABS(HYD(NZ,NDOLD)) .LT. 0.1) THEN
	  CONST=1.0D0
	ELSE
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Possible error on REGRID_B_ON_NE - b(n) not equal'
	  WRITE(LUER,*)'to 1.0/pm0.1 or 0.0/pm0.1 at depth'
	END IF
	CONST=0.0D0
!
! Determine radius at which optical depth is unity. We use NEW_ED as a
! temporary vector for ED*CLUMP_FAC.
!
        NEW_ED(1:ND)=ED(1:ND)*CLUMP_FAC(1:ND)
        TAU(1)=6.65D-15*NEW_ED(1)*R(1)
        K=1
	DO I=2,ND
          TAU(I)=TAU(I-1)+6.65D-15*(NEW_ED(I-1)+NEW_ED(I))*(R(I-1)-R(I))*0.5D0
          IF(TAU(I) .LE. 1.0D0)K=I
        END DO
	IF(K .EQ. 1 .OR. K .EQ. ND)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error computing RTAU1 in REGRID_B_ON_SPH_TAU'
        END IF
	T1=(1.0D0-TAU(K))/(TAU(K+1)-TAU(K))
        RTAU1=T1*R(K+1)+(1.0D0-T1)*R(K)
!
! Compute the spherical optical depth scale. Assume atmosphere has constant
! V at outer boundary.
!
        DO I=1,ND
          NEW_ED(I)=NEW_ED(I)*RTAU1*RTAU1/R(I)/R(I)
        END DO
        TAU(1)=6.65D-15*NEW_ED(1)*R(1)/3.0D0
        DO I=2,ND
          TAU(I)=TAU(I-1)+6.65D-15*(NEW_ED(I-1)+NEW_ED(I))*(R(I-1)-R(I))*0.5D0
        END DO
!
! Determine radius at which optical depth is unity in old model. OLDED has already
! been modified by the clumping factor.
!
        OLDTAU(1)=6.65D-15*OLDED(1)*ROLD(1)
        K=1
	DO I=2,NDOLD
          OLDTAU(I)=OLDTAU(I-1)+6.65D-15*(OLDED(I-1)+OLDED(I))*(ROLD(I-1)-ROLD(I))*0.5D0
          IF(OLDTAU(I) .LE. 1.0D0)K=I
        END DO
	IF(K .EQ. 1 .OR. K .EQ. ND)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error computing RTAU1_OLD in REGRID_B_ON_SPH_TAU'
        END IF
        T1=(1.0D0-OLDTAU(K))/(OLDTAU(K+1)-OLDTAU(K))
        RTAU1_OLD=T1*ROLD(K+1)+(1.0D0-T1)*ROLD(K)
!
! Compute the spherical optical depth scale.
!
        DO I=1,NDOLD
          OLDED(I)=OLDED(I)*(RTAU1_OLD/ROLD(I))**2
        END DO
        OLDTAU(1)=6.65D-15*OLDED(1)*ROLD(1)/3.0D0
        DO I=2,NDOLD
          OLDTAU(I)=OLDTAU(I-1)+6.65D-15*(OLDED(I-1)+OLDED(I))*(ROLD(I-1)-ROLD(I))*0.5D0
        END DO
!
! Detrmine whether new mesh extends beyond oldmesh.
! NX is used to define the region over which we may use a linear
! interpolation. Beyond the old radius mesh we set the departure coefficients
! constant or equal to 1.
!
	NXST=1
	DO WHILE (TAU(NXST) .LT. OLDTAU(1))
	  NXST=NXST+1
	END DO
	NX_END=ND
	DO WHILE (TAU(NX_END) .GT. OLDTAU(NDOLD))
	  NX_END=NX_END-1
	END DO
	NX=NX_END-NXST+1
!
! Interpolate the populations assuming that the departure
! coefficients remain constant. Constant is one if read in b's-1.0
!
	DO J=1,NZ
	  DO I=1,NDOLD
	    TB(I)=LOG(HYD(J,I)+CONST)
	  END DO
	  CALL LINPOP(TAU(NXST),TA(NXST),NX,OLDTAU,TB,NDOLD)
	  DO I=NXST,NX_END
	    DHEN(J,I)=EXP(TA(I))
	  END DO
!
! Assume d.c. is constant in outer region. This option is consistent with 
! the assumption of constant T.
!
	  IF(NXST .NE. 1)THEN
	    DO I=1,NXST-1
	      DHEN(J,I)=EXP(TB(1))
	    END DO
	  END IF
!
	  IF(NX_END .LT. ND)THEN
	    DO I=ND,NX_END+1,-1
	      DHEN(J,I)=DHEN(J,NX_END)
	    END DO
	  END IF
	END DO
!
! Compute departure coefficints for N>NZ . New levels are set to have the same
! departure coefficient as highest known level.
!
	IF(N .GT. NZ)THEN
	  DO I=1,ND
	    DO J=NZ+1,N
	      DHEN(J,I)=DHEN(NZ,I)
	    END DO
	  END DO
	END IF
!
! Regrid ion density.
!
	CALL LINPOP(TAU(NXST),ION(NXST),NX,OLDTAU,OLDION,NDOLD)
	DO I=1,NX_END
	  ION(I)=EXP(ION(I))
	END DO
!
! Assume ion density is fixed relative to ion density.
!
	IF(NX_END .LT. ND)THEN
	  DO I=ND,NX_END+1,-1
	    ION(I)=ED(I)*ION(NX_END)/ED(NX_END)
	  END DO
	END IF
!
! Free memory.
!
	DEALLOCATE (HYD)
	DEALLOCATE (ROLD)
	DEALLOCATE (OLDED)
	DEALLOCATE (OLDTAU)
	DEALLOCATE (TA)
	DEALLOCATE (TB)
	DEALLOCATE (OLDION)
!
	RETURN
	END
