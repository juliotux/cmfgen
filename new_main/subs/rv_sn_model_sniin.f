!
!   14 July 2014
!   Routine to set up the radius grid for IIn simulations. We use the TAU_ES
!   grid from hydrodynamical simulation
!
	SUBROUTINE RV_SN_MODEL_SNIIN(R,V,SIGMA,RMAX,RP,VCORE,BETA1,RDINR,LU,ND)
	IMPLICIT NONE
	INTEGER ND
	INTEGER LU
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 SIGMA(ND)
!
	REAL*8 RMAX
	REAL*8 RP
	REAL*8 VCORE
	REAL*8 BETA1
	LOGICAL RDINR
!
! Local arrays.
!
	REAL*8 TA(ND),TB(ND),TC(ND)
	REAL*8 T1,T2,DLNR
!
	INTEGER NBND_INS
	INTEGER I,J,MND
	INTEGER IOS,NOLD,NDOLD
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	CHARACTER*80 STRING
!
         INTEGER NHYDRO
         REAL*8 XV,XD,XT,XL,XTAU,RAT
         REAL*8 NEW_TAU(ND)
         REAL*8, ALLOCATABLE :: R_HYDRO(:),TAU_HYDRO(:),V_HYDRO(:)

! Variables for acceleration zone
         LOGICAL ADD_ACC_ZONE
         REAL*8 T0, S1, S2, RINT, HRHO, VRAT, VMIN, BETA_ACC
!
!
       NBND_INS=3 ! LUC: changed from 2 to 3 --- KEYWORD in VADAT not used
!
! Check whether the passed parameters are valid.
!
	IF(BETA1 .LT. 0.0D0)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error in RV_SN_MODEL --- Invalid BETA'
	  STOP
	END IF
	IF(NBND_INS .LT. 1 .OR. NBND_INS .GT. 3)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error in STARPCYG_V3 --- Invalid NBND_INS'
          WRITE(LUER,*)'NBND_INS should be 1, 2, or 3'
	END IF
!
	IF(RDINR)THEN
	  OPEN(UNIT=LU,STATUS='OLD',FILE='RDINR',IOSTAT=IOS)
          IF(IOS .NE. 0)THEN
            LUER=ERROR_LU()
            WRITE(LUER,*)'Error in RV_SN_MODEL_02 --- File with R grid not found'
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
	    WRITE(LUER,*)'Error in RV_SN_MODEL_02 --- unable to read header in file with R grid'
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
	      WRITE(LUER,*)'Error in RV_SN_MODEL_02 --- unable to read R grid from file'
	      STOP
	    END IF
	  END DO
	  R(1)=RMAX
!
! Compute Velocity and SIGMA
!
	  DO I=1,ND
	    V(I)=VCORE*(R(I)/R(ND))**BETA1
	    SIGMA(I)=BETA1-1.0D0
	  END DO
	  R(ND)=RP
	  CLOSE(UNIT=LU)
	  RETURN
	END IF
!
!   Read the hydro input file
!
        OPEN(UNIT=LU,STATUS='OLD',FILE='input_hydro.dat',IOSTAT=IOS)
        IF(IOS .NE. 0)THEN
            LUER=ERROR_LU()
            WRITE(LUER,*)'Error in RV_SN_MODEL_SNIIN --- File input_hydro.dat not found'
            STOP
        END IF
        READ(LU,*) ADD_ACC_ZONE
        IF (ADD_ACC_ZONE) THEN
             READ(LU,*) RINT
             READ(LU,*) HRHO
             READ(LU,*) VRAT
             READ(LU,*) VMIN
             READ(LU,*) BETA_ACC
        ENDIF
        READ(LU,*) STRING
        READ(LU,*) NHYDRO
        READ(LU,*) STRING  ! reads the bogus line from python script
! TB is the radius and TC is tau_es scaled with the temperature to resolve shock
        ALLOCATE (R_HYDRO(NHYDRO),V_HYDRO(NHYDRO),TAU_HYDRO(NHYDRO))
        DO I=1,NHYDRO
           READ(LU,*) R_HYDRO(I),V_HYDRO,XD,XT,XL,XTAU,TAU_HYDRO(I)
        ENDDO
        CLOSE(LU)

        MND=ND-2*NBND_INS
        T1 = TAU_HYDRO(1)
        T2 = TAU_HYDRO(NHYDRO)
        RAT = EXP( LOG(T2/T1)/(MND) )

        DO I=1,MND
          NEW_TAU(I) = T1 * RAT**(I-1)
       ENDDO
       NEW_TAU(1) = TAU_HYDRO(1)
       NEW_TAU(MND) = TAU_HYDRO(NHYDRO)
       DO I=1,NHYDRO
           WRITE(127,*) R_HYDRO(I),TAU_HYDRO(I)
       ENDDO
       DO I=1,MND
           WRITE(127,*) NEW_TAU(I)
       ENDDO
       CALL FLUSH(127)
!
! Change TAU_HYDRO and NEW_TAU to a log to improve the interpolation.
!
       TAU_HYDRO = LOG(TAU_HYDRO)
       NEW_TAU   = LOG(NEW_TAU)

        CALL LIN_INTERP(NEW_TAU,TA,MND,TAU_HYDRO,R_HYDRO,NHYDRO)

        ! reverse order of the TA array to become hte R grid
        TB = TA
        DO I=1,MND
           TA(I) = TB(MND-I+1)
        ENDDO

        RP = TA(MND)
        RMAX=  TA(1)
        WRITE(127,*) 'LUC: NBND_INS: ',NBND_INS
        DO I=1,MND
           WRITE(127,'(ES14.4)')TA(I)
        ENDDO
                CALL FLUSH(127)
!
!	MND=ND-2*NBND_INS
!	T1=LOG(RMAX/RP)
!	T1=EXP(T1/(MND-1))
!	TA(MND)=RP
!	DO I=MND-1,2,-1
!	  TA(I)=RP*(T1**(MND-I))
!	  WRITE(127,'(ES14.4)')TA(I)
!	END DO
!	TA(1)=RMAX
!
! Insert finer grid near both boundaries.
!
	DO I=2,MND-1
	  R(I+NBND_INS)=TA(I)
	END DO
	R(1)=TA(1)
	R(ND)=TA(MND)
	IF(NBND_INS .EQ. 1)THEN
	  R(2)=TA(1)-(TA(1)-TA(2))/20.0
	  R(ND-1)=R(ND)+(TA(MND-1)-TA(MND))/20.0D0
	ELSE IF(NBND_INS .EQ. 2)THEN
!	  R(2)=TA(1)-(TA(1)-TA(2))/10.0D0
!	  R(ND-1)=R(ND)+(TA(MND-1)-TA(MND))/10.0D0
!	  R(3)=TA(1)-(TA(1)-TA(2))/3.0D0
!	  R(ND-2)=R(ND)+(TA(MND-1)-TA(MND))/3.0D0
	  R(2)=TA(1)-(TA(1)-TA(2))/50.0D0
	  R(ND-1)=R(ND)+(TA(MND-1)-TA(MND))/10.0D0
	  R(3)=TA(1)-(TA(1)-TA(2))/5.0D0
	  R(ND-2)=R(ND)+(TA(MND-1)-TA(MND))/3.0D0
	ELSE IF(NBND_INS .EQ. 3)THEN
	  R(2)=TA(1)-(TA(1)-TA(2))/20.0D0
	  R(ND-1)=R(ND)+(TA(MND-1)-TA(MND))/20.0D0
	  R(3)=TA(1)-(TA(1)-TA(2))/8.0D0
	  R(ND-2)=R(ND)+(TA(MND-1)-TA(MND))/8.0D0
	  R(4)=TA(1)-(TA(1)-TA(2))/3.0D0
	  R(ND-3)=R(ND)+(TA(MND-1)-TA(MND))/3.0D0
	END IF
!
	DO I=1,NHYDRO
	  T1=R_HYDRO(I)
	  R_HYDRO(I)=R_HYDRO(NHYDRO-I+1)
	  R_HYDRO(NHYDRO-I+1)=T1
	  T1=V_HYDRO(I)
	  V_HYDRO(I)=V_HYDRO(NHYDRO-I+1)
	  V_HYDRO(NHYDRO-I+1)=T1
	END DO
!
        CALL LIN_INTERP(V,TA,ND,V_HYDRO,R_HYDRO,NHYDRO)
	WRITE(127,*)RP,RMAX,VCORE,BETA1
	DO I=2,ND-1
	  T1=(V(I-1)-V(I))/(R(I-1)-R(I))
	  T2=(V(I)-V(I+1))/(R(I)-R(I+1))
	  SIGMA(I)=(T1*(R(I)-R(I+1))+T2*(R(I-1)-R(I)))/(R(I-1)-R(I))
	  SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0D0
	END DO
	SIGMA(1)=R(1)*(V(1)-V(2))/(R(1)-R(2))/V(1)-1.0D0
	SIGMA(ND)=R(ND)*(V(ND-1)-V(ND))/(R(ND-1)-R(ND))/V(ND)-1.0D0
!
	DO I=1,ND
	  WRITE(127,'(3ES12.4)')R(I),V(I),SIGMA(I)
	END DO
        CALL FLUSH(127)
        DEALLOCATE (R_HYDRO,V_HYDRO,TAU_HYDRO)
	RETURN
!
!
! Compute Velocity and SIGMA
!
       IF (ADD_ACC_ZONE) THEN
         HRHO = HRHO * R(ND)
         DO I=1,ND
           IF (R(I).GT.RINT) THEN
               V(I)=VCORE*(R(I)/RINT)**BETA1
               SIGMA(I)=BETA1-1.0D0
           ELSE
               T0 = (VCORE-VMIN) / (1.-R(ND)/RINT)**BETA_ACC ! V=VCORE at RINT
               IF (I.EQ.ND) THEN
                   T1 = VMIN
               ELSE
                   T1 = VMIN + T0   * (1.-R(ND)/R(I))**BETA_ACC
               ENDIF
               T2 = 1.    + VRAT * EXP( (R(ND)-R(I))/HRHO )
               V(I) = T1 /  T2


               IF (I.EQ.ND) THEN
                   S1 = 0.
               ELSE
                   S1 = T0*BETA_ACC*R(ND)*(1.-R(ND)/R(I))**(BETA_ACC-1.)/R(I)/T1
               ENDIF
               S2 = R(I)*VRAT*EXP( (R(ND)-R(I))/HRHO ) / T2 / HRHO
               SIGMA(I) = S1 + S2 - 1.0

           ENDIF
        ENDDO
        ! Patch to avoid rapid change in SIGMA at base
        !SIGMA(ND) = SIGMA(ND-1)
       ELSE
          DO I=1,ND
	        V(I)=VCORE*(R(I)/R(ND))**BETA1
	        SIGMA(I)=BETA1-1.0D0
	      END DO
       ENDIF
!
	WRITE(127,*)RP,RMAX,VCORE,BETA1
	DO I=1,ND
	  WRITE(127,'(3ES12.4)')R(I),V(I),SIGMA(I)
	END DO
         CALL FLUSH(127)
!
	RETURN
	END
