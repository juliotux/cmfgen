!
! Subroutine to read a NYDRO model from SN_HYDRO_DATA and genrate
! the R grid, and the V and SIGMA vectors.
!
	SUBROUTINE SET_RV_HYDRO_MODEL(R,V,SIGMA,SN_AGE_DAYS,RMAX,RCORE,RDINR,ND,LU)
	IMPLICIT NONE
!
! Created 14-Sep-2006
!
	INTEGER ND
	INTEGER LU
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 SIGMA(ND)
	REAL*8 SN_AGE_DAYS
	REAL*8 RMAX
	REAL*8 RCORE
	LOGICAL RDINR
!
! Local arrays and variables.
!
	REAL*8, ALLOCATABLE :: R_HYDRO(:)
	REAL*8, ALLOCATABLE :: V_HYDRO(:)
	REAL*8, ALLOCATABLE :: SIGMA_HYDRO(:)
	REAL*8, ALLOCATABLE :: DENSITY_HYDRO(:)
	REAL*8, ALLOCATABLE :: KAPPA_HYDRO(:)
	REAL*8, ALLOCATABLE :: TAU_HYDRO(:)
	REAL*8, ALLOCATABLE :: LOG_R_HYDRO(:)
!
	REAL*8 TA(ND)
	REAL*8 XN(ND)
	REAL*8 TAU(2*ND)
	REAL*8 LOG_R(2*ND)
	REAL*8 OLD_SN_AGE_DAYS
	REAL*8 SN_EXP_FACTOR
	REAL*8 dTAU
	REAL*8 dLOGR
	REAL*8 T1
	REAL*8 NEXT_R
	REAL*8 TAU_BEG,TAU_END
	INTEGER ND_TMP
	INTEGER NX
	INTEGER NSP
	INTEGER IOS
	INTEGER NOLD,NDOLD
	INTEGER I,J,L
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	INTEGER, PARAMETER :: IONE=1
	CHARACTER(LEN=200) STRING
!
	WRITE(6,*)'Entering SET_RV_HYDRO',RDINR
	OPEN(UNIT=LU,FILE='SN_HYDRO_DATA',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in SET_RV_HYDRO_MODEL --- SN_HYDRO_DATA file not found'
	    WRITE(LUER,*)'Create file or EDIT option in VADAT'
	    STOP
	  END IF
!
! Get the number of data points in the HYDRO model, and the number 
! of species.
!
	NX=0; NSP=0; OLD_SN_AGE_DAYS=0.0D0
	DO WHILE (1 .EQ. 1)
	  STRING=' '
	  DO WHILE (STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	    READ(LU,'(A)')STRING 
	  END DO
	  IF(INDEX(STRING,'Number of data points:') .NE. 0)THEN
	    I=INDEX(STRING,'points:')+7
	    READ(STRING(I:),*)NX
	  ELSE IF(INDEX(STRING,'Number of mass fractions:') .NE. 0)THEN
	    I=INDEX(STRING,'fractions:')+10
	    READ(STRING(I:),*)NSP
	  ELSE IF(INDEX(STRING,'Time(days) since explosion:') .NE. 0)THEN
	    I=INDEX(STRING,'explosion:')+10
	    READ(STRING(I:),*)OLD_SN_AGE_DAYS
	  ELSE IF(INDEX(STRING,'Radius grid') .NE. 0)THEN
	     IF(NSP .EQ. 0 .OR. NX .EQ. 0 .OR. OLD_SN_AGE_DAYS .EQ. 0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in RD_SN_DATA'
	      WRITE(LUER,*)'NSP or NX undefined'
	      STOP
	    ELSE
	      EXIT
	    END IF
	  END IF
	END DO
	WRITE(6,*)NX,NSP,OLD_SN_AGE_DAYS
!
	SN_EXP_FACTOR=SN_AGE_DAYS/OLD_SN_AGE_DAYS
!
	ALLOCATE (V_HYDRO(NX));          V_HYDRO=0.0D0
	ALLOCATE (SIGMA_HYDRO(NX));      SIGMA_HYDRO=0.0D0
	ALLOCATE (DENSITY_HYDRO(NX));    DENSITY_HYDRO=0.0D0
	ALLOCATE (KAPPA_HYDRO(NX));      KAPPA_HYDRO=0.0D0
	ALLOCATE (R_HYDRO(NX));          R_HYDRO=0.0D0
	ALLOCATE (TAU_HYDRO(NX));        TAU_HYDRO=0.0D0
	ALLOCATE (LOG_R_HYDRO(NX));      LOG_R_HYDRO=0.0D0
!
! Get basic HYDRO grid vectors.

	DO WHILE (1 .EQ. 1)
	  DO WHILE (STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	    READ(LU,'(A)')STRING 
	  END DO
	  IF(INDEX(STRING,'Radius grid') .NE. 0)THEN
	    READ(LU,*)R_HYDRO
	    R_HYDRO=R_HYDRO*SN_EXP_FACTOR
	  ELSE IF(INDEX(STRING,'Velocity') .NE. 0)THEN
	    READ(LU,*)V_HYDRO
	  ELSE IF(INDEX(STRING,'Sigma') .NE. 0)THEN
	    READ(LU,*)SIGMA_HYDRO
	  ELSE IF(INDEX(STRING,'Density') .NE. 0)THEN
	    READ(LU,*)DENSITY_HYDRO
	    DENSITY_HYDRO=DENSITY_HYDRO/(SN_EXP_FACTOR**3)
	  ELSE IF(INDEX(STRING,'Kappa') .NE. 0)THEN
	    READ(LU,*)KAPPA_HYDRO
	  ELSE IF(INDEX(STRING,'ass fraction') .NE. 0)THEN
	    IF(R_HYDRO(1) .EQ. 0.0D0 .OR. 
	1             V_HYDRO(1) .EQ. 0.0D0 .OR.
	1             DENSITY_HYDRO(1) .EQ. 0.0D0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error reading SN data in SET_RV_HYDRO_MODEL'
	      WRITE(LUER,*)'R, V, the DENSITY is zero'
	      STOP
	    ELSE IF(KAPPA_HYDRO(1) .EQ. 0.0D0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Warning: reading SN data in SET_RV_HYDRO_MODEL'
	      WRITE(LUER,*)'KAPPA is zero'
	      EXIT
	    ELSE
	      EXIT
	    END IF
	  END IF
	  STRING=' '
	END DO
	WRITE(6,*)'Successfuly read data in SET_RV_HYDRO'
	CLOSE(LU)
!
	IF(RDINR)THEN
	  OPEN(UNIT=LU,STATUS='OLD',FILE='RDINR',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in SET_RV_HYDRO_MODEL --- File with R grid not found'
	    WRITE(LUER,*)'Create file or EDIT option in VADAT'
	    STOP
	  END IF
!
! Check whether the file has a record containing 'Format date'. Its presence
! effects the way we read the file.
!
	  I=0
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	    I=I+1
	    READ(LU,'(A)')STRING
	  END DO
	  IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(LU)
!
	  READ(LU,*,IOSTAT=IOS)TA(1),TA(1),NOLD,NDOLD
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in SET_RV_HYDRO_MODEL --- unable to read header in file with R grid'
	    STOP
	  END IF
!
! Check relative values.
!
	  IF(ND .NE. NDOLD)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error-NDOLD and ND are not equal in RDINR'
	    WRITE(LUER,*)'NDOLD=',NDOLD,' ND=',ND
	    STOP
	  END IF
!
! TA is used for everything but R which is all we want.
!
	  DO I=1,ND
	    READ(LU,*,IOSTAT=IOS)R(I),TA(I),TA(I),TA(I)
	    IF(IOS .EQ. 0)READ(LU,*,IOSTAT=IOS)(TA(J),J=1,NOLD)
	    IF(IOS .NE. 0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in SET_RV_HYDRO_MODEL --- unable to R grid in file with R grid'
	      STOP
	    END IF
          END DO
	  DO I=1,ND               	!need to issue a warning here
	    R(I)=R_HYDRO(NX)*(R(I)/R(ND))
	  END DO
	  RMAX=R(1); RCORE=R(ND)
	  R(1)=MIN(R(1),R_HYDRO(1))
!
! Since Hubble law.
!
	  V=R*(V_HYDRO(1)/R(1))
	  SIGMA=0.0D0
	  WRITE(6,*)'Read in R grid from RDINR'
!
	ELSE
!
! Compute TAU. Since we are only setting the grid, a simple trapzoidal
! rule integration will suffice.
!
	  KAPPA_HYDRO=1.0D+10*KAPPA_HYDRO*DENSITY_HYDRO
	  TAU_HYDRO(1)=KAPPA_HYDRO(1)*R_HYDRO(1)
	  T1=LOG(KAPPA_HYDRO(2)/KAPPA_HYDRO(1))/LOG(R_HYDRO(1)/R_HYDRO(2))
	  IF(T1 .GT. 1.5D0)TAU_HYDRO(1)=TAU_HYDRO(1)/(T1-1)
	  DO I=2,NX
	    TAU_HYDRO(I)=TAU_HYDRO(I-1)+0.5D0*(KAPPA_HYDRO(I-1)+KAPPA_HYDRO(I))*
	1                  (R_HYDRO(I-1)-R_HYDRO(I))
	  END DO
!
! Estimate spacing to get required grid spacing.
!
	  TAU_HYDRO=LOG(TAU_HYDRO)
	  dTAU=(TAU_HYDRO(NX)-TAU_HYDRO(1))/(ND-7)
	  LOG_R_HYDRO=LOG(R_HYDRO)
	  dLOGR=(LOG_R_HYDRO(1)-LOG_R_HYDRO(NX))/(ND-7)
	  DO I=1,NX
	    WRITE(177,*)I,LOG_R_HYDRO(I),TAU_HYDRO(I)
	  END DO
	  WRITE(177,*)'dTAU=',dTAU
	  WRITE(177,*)'dLOGR=',dLOGR
!
	  LOG_R(1)=LOG_R_HYDRO(1)
	  LOG_R(2)=LOG_R_HYDRO(1)-0.01*dLOGR
	  LOG_R(3)=LOG_R_HYDRO(1)-0.03*dLOGR
	  LOG_R(4)=LOG_R_HYDRO(1)-0.1*dLOGR
	  LOG_R(5)=LOG_R_HYDRO(1)-0.4*dLOGR
	  LOG_R(6)=LOG_R_HYDRO(1)-dLOGR
	  J=1; I=6
	  DO WHILE(1. EQ. 1)
	    I=I+1
	    IF(I .GT. 2*ND-4)THEN
	      WRITE(LUER,*)'Error in SET_RV_HYDRO_MODEL --- LOG_R and TAU vectors too small'
	      STOP
	    END IF
	    DO WHILE(LOG_R_HYDRO(J+1) .GT. LOG_R(I-1))
	      J=J+1
	    END DO
	    T1=(LOG_R_HYDRO(J)-LOG_R(I-1))/(LOG_R_HYDRO(J)-LOG_R_HYDRO(J+1))
	    TAU_BEG=T1*TAU_HYDRO(J+1)+(1.0D0-T1)*TAU_HYDRO(J)
	    NEXT_R=LOG_R(I-1)-dLOGR
	    DO WHILE(LOG_R_HYDRO(J+1) .GT. NEXT_R)
	      J=J+1
	    END DO
	    T1=(LOG_R_HYDRO(J)-NEXT_R)/(LOG_R_HYDRO(J)-LOG_R_HYDRO(J+1))
	    TAU_END=T1*TAU_HYDRO(J+1)+(1.0D0-T1)*TAU_HYDRO(J)
	    IF(TAU_END-TAU_BEG .GT. dTAU)THEN
	      NEXT_R=LOG_R(I-1)-dLOGR*dTAU/(TAU_END-TAU_BEG)
	    END IF
	    LOG_R(I)=NEXT_R
	    IF(LOG_R(I)-1.5*dLOGR .LT. LOG_R_HYDRO(NX))EXIT
	    WRITE(177,*)I,LOG_R(I)
	  END DO
	  dLOGR=LOG_R(I)-LOG_R_HYDRO(NX)
	  LOG_R(I+3)=LOG_R_HYDRO(NX)
	  LOG_R(I+2)=LOG_R_HYDRO(NX)+0.1*dLOGR
	  LOG_R(I+1)=LOG_R_HYDRO(NX)+0.4*dLOGR
	  ND_TMP=I+3
!
	  DO I=1,ND_TMP; TAU(I)=I; END DO
	  DO I=1,ND
	    XN(I)=1.0D0+(I-1.0D0)*(ND_TMP-1.0D0)/(ND-1.0D0)
	  END DO
	  DO I=1,ND_TMP
	    WRITE(178,*)I,LOG_R(I),TAU(I)
	  END DO
	  CALL MON_INTERP(R,ND,IONE,XN,ND,LOG_R,ND_TMP,TAU,ND_TMP)
	  DO I=1,ND
	    WRITE(177,*)I,XN(I),R(I)
	  END DO
!
	  R=EXP(R)
	  R(1)=R_HYDRO(1)
	  R(ND)=R_HYDRO(NX)
	  CLOSE(UNIT=177)
!
! Since Hubble law.
!
	  V=R*(V_HYDRO(1)/R_HYDRO(1))
	  SIGMA=0.0D0
	  RMAX=R(1); RCORE=R(ND)
	  WRITE(6,*)'Computed R grid'
!
	END IF
!
	DEALLOCATE (R_HYDRO, V_HYDRO, SIGMA_HYDRO)
	DEALLOCATE (DENSITY_HYDRO, KAPPA_HYDRO, TAU_HYDRO, LOG_R_HYDRO)
	WRITE(6,*)'Exiting SET_RV_HYDRO'
!
	RETURN
	END
