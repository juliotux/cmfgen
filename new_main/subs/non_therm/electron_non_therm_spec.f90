	SUBROUTINE ELECTRON_NON_THERM_SPEC(ND)
	USE MOD_NON_THERM
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
! Altered 15-Jun-2016 : NT_ITERATION_COUNTER is now defined in CONTROL_VARIABLE_MOD
! Altered 23-Feb-2102 : ASUM was changed to QSUM, since A can be zero for some transition.
! Altered 16-Feb-2012 : Cleaned (many prior changes to this date also).
! Altered  1-Dec-2011 : Included NON_THERM_IT_CNTRL and NT_OMIT_ION_SCALE.
!                       NT_LEV_SCALE changed to NT_OMIT_LEV_SCALE
!                       Altered method for excluding levels.
!                       Fixed bugs with excitaton sction.
! Altered  7-Nov-2011 : Altered selection of levels used for evaluating non-thermal spectrum.
!                       We now use the atom density, rather than ion density at tthe depth of interest.
!                       NT_OMIT_LEV_SCALE added 23-Nov-2011
!                       CONTROL_VARIABL_MOD included.
!                       Set SOURCE to 0 before setting its value: important for CONSTANT and BELL_SHAPE
!                       Installed NT_SOURCE_TYPE variable.
!
!	INTEGER NUM_IONS
	INTEGER ND
!
	REAL*8, SAVE, ALLOCATABLE :: LELEC(:)
	REAL*8, SAVE, ALLOCATABLE :: SOURCE(:)
	REAL*8, SAVE, ALLOCATABLE :: RHS(:)
	REAL*8, SAVE, ALLOCATABLE :: DETAI_DE(:)
	REAL*8, SAVE, ALLOCATABLE :: MAT(:,:)
	REAL*8, SAVE, ALLOCATABLE :: Qnn(:)
	REAL*8, SAVE, ALLOCATABLE :: Qnn_TMP(:)
	REAL*8, SAVE, ALLOCATABLE :: POP_VEC(:)
!
	INTEGER, SAVE, ALLOCATABLE :: INDX(:)
	LOGICAL, SAVE, ALLOCATABLE :: DO_THIS_ION(:)
	logical, save, allocatable :: do_and_1st_time(:)
	integer, save, allocatable :: nlow_maxs(:)
!
	REAL*8 ION_SUM(NUM_IONS)
	REAL*8 SPEC_SUM(NUM_IONS)
	LOGICAL DO_THIS_ION_EXC(NUM_IONS)
!
	REAL*8 EKT
	REAL*8 EKTP
	REAL*8 EMIN
	REAL*8 EMAX
	REAL*8 E_INIT
	REAL*8 dE
!
	REAL*8 T1,T2
	REAL*8 RATE
	REAL*8 SIG1,SIG2
	REAL*8 XCROSS
	REAL*8 NA_dX_CROSS
	REAL*8 XION_POT
	REAL*8 NATOM
	REAL*8 QSUM
	REAL*8 SCALER_SPEC_SUM
	REAL*8 SCALER_ION_SUM
	REAL*8 dXKT_MIN
!
	INTEGER GET_INDX_DP
	REAL*8 INTSIGC
	EXTERNAL INTSIGC,GET_INDX_DP
!
! We now set these values by values passed from CMFGEN.
!
	REAL*8 XKT_MIN
	REAL*8 XKT_MAX
!
! These parameters (except Hz_TO_eV) probably only need to be changed if performing tests.
!
	REAL*8, PARAMETER :: Hz_to_eV=13.60569253D0/3.289841960D0
	REAL*8, PARAMETER :: DELTA_ENR_SOURCE=30.0D0
	LOGICAL, PARAMETER :: INCLUDE_EXCITATION=.TRUE.
	LOGICAL, PARAMETER :: INCLUDE_IONIZATION=.TRUE.
	CHARACTER(LEN=3), PARAMETER :: XKT_METHOD='lin'
!
	INTEGER ISPEC
	INTEGER ID
	INTEGER IT
	INTEGER IKT
	INTEGER IKT0
	INTEGER IKTP
	INTEGER IKTPN
	INTEGER DPTH_INDX
	INTEGER I, J, K
	INTEGER IST
	INTEGER NL,NUP
	INTEGER L
	INTEGER MAX_LOW_LEV
!
	INTEGER, PARAMETER :: LU_ER=6
	INTEGER, PARAMETER :: LU_TH=8
	integer, parameter :: LU_BETHE = 100
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
	LOGICAL DO_LOC_CHK
!
	LOGICAL INJECT_DIRAC_AT_EMAX
	CHARACTER(LEN=10) SOURCE_TYPE
	CHARACTER(LEN=40) TMP_NAME
!
!
	XKT_MIN=NT_EMIN	
	XKT_MAX=NT_EMAX
	INJECT_DIRAC_AT_EMAX=.FALSE.
	IF(NT_SOURCE_TYPE(1:12) .EQ. 'INJECT_DIRAC')THEN
	  INJECT_DIRAC_AT_EMAX=.TRUE.
	ELSE IF (NT_SOURCE_TYPE .EQ. 'CONSTANT')THEN
	  SOURCE_TYPE='CONSTANT'
	ELSE IF (NT_SOURCE_TYPE .EQ. 'BELL_SHAPE')THEN
	  SOURCE_TYPE='BELL_SHAPE'
	ELSE
	  WRITE(LU_ER,*)'Error in electron_non_therm_spec.f90 -- SOURCE_TYPE not recognized'
	  WRITE(LU_ER,*)'SOURCE_TYPE= ',TRIM(SOURCE_TYPE)
	  STOP 
	END IF
!
	IF(FIRST_TIME)THEN
	  NT_ITERATION_COUNTER=1
	ELSE
	  NT_ITERATION_COUNTER=NT_ITERATION_COUNTER+1
	  IF(MOD(NT_ITERATION_COUNTER-1,NON_THERMAL_IT_CNTRL) .NE. 0)RETURN
	END IF
	DO_LOC_CHK=.TRUE.
!
	WRITE(LU_ER,*)' '
	WRITE(LU_ER,*)'Entering ELECTRON_NON_THERM_SPEC'
	OPEN(UNIT=LU_TH,FILE='NON_THERM_SPEC_INFO',STATUS='UNKNOWN')
	CALL SET_LINE_BUFFERING(LU_TH)
	WRITE(LU_TH,'(A)')' '
	WRITE(LU_TH,'(A,I5)')  ' Number of enrgy bins is ',NKT
	WRITE(LU_TH,'(A,F8.2)')' Emin in eV is',XKT_MIN
	WRITE(LU_TH,'(A,F8.2)')' Emax in eV is',XKT_MAX
	WRITE(LU_TH,'(A)')' '
!
! Allocate arrays that will be used for each depth:.
!
! The following vectors/arrays are stored in MOD_NON_THERM
!
	IF(.NOT. ALLOCATED(XKT))THEN
	  NKT=NT_NKT
	  WRITE(LU_TH,*)'Allocating memory to describe non-thermal electron distribution'
	  ALLOCATE (XKT(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (dXKT(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (dXKT_ON_XKT(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (YE(NKT,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (FRAC_ELEC_HEATING(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (FRAC_ION_HEATING(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (FRAC_EXCITE_HEATING(ND),STAT=IOS)
	  IF (IOS.NE.0) THEN
            WRITE(LU_ER,*) 'Error allocating XKT',IOS
	    STOP
	  END IF
	  XKT(1:NKT)=0.0D0; dXKT(1:NKT)=0.0D0
	END IF
!
! The following vectors/arrrays are local:
!
! LELEC : Energy loss through Coulomb interaction
! SOURCE : Electron source
!
	WRITE(LU_TH,*)'Checking MAT memory allocation'
	IF(.NOT. ALLOCATED(MAT))THEN
	  WRITE(LU_TH,*)'Allocating memory in ELECTRON_NON_THERMAL_SPEC'
	  ALLOCATE (LELEC(NKT),STAT=IOS) 
	  IF(IOS .EQ. 0)ALLOCATE (SOURCE(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RHS(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (Qnn(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (Qnn_TMP(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (dETAI_dE(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (MAT(NKT,NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (POP_VEC(NUM_THD))
	  IF(IOS .EQ. 0)ALLOCATE (INDX(NUM_THD))
	  IF(IOS .EQ. 0)ALLOCATE (DO_THIS_ION(NUM_IONS))
	  if(ios .eq. 0)allocate (do_and_1st_time(num_ions))
	  if(ios .eq. 0)allocate (nlow_maxs(num_ions))
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_ER,*) 'Error allocating LELEC, SOURCE, etc in ELECTRON_NON_TERMAL_SPEC',IOS
	    STOP
	  END IF
	END IF
!
! Set the energy array. Use a linear scaling since we need to resolve the high-energy source
! but also the transition near the ionization/excitation potentials at about 10eV.
!
	IF(FIRST_TIME)THEN
!
	  WRITE(LU_TH,*)'Setting XKT array'
	  CALL SET_XKT_ARRAY(XKT_MIN,XKT_MAX,NKT,XKT,dXKT,XKT_METHOD)
	  dXKT_ON_XKT(1:NKT)=dXKT(1:NKT)/XKT(1:NKT)
!
! If needed, allocate memory for collisional ioinzation cross-sections.
!
	  IOS=0
	  DO IT=1,NUM_THD
	    IF(THD(IT)%PRES)ALLOCATE (THD(IT)%CROSS_SEC(NKT),STAT=IOS) 
	    IF(IOS.NE.0) THEN
              WRITE(LU_ER,*) 'Error allocating THD%CROSS_SEC',IOS
	      STOP
	    END IF
	  END DO
	  WRITE(LU_TH,*)'Computing cross-sections'
	  CALL ARNAUD_CROSS_V3()
!
	  I=0
	  DO IT=1,NUM_THD
	    IF(THD(IT)%PRES)I=I+1
	  END DO
	  WRITE(LU_TH,*)'Number of active species in ELECTRON_NON_TERMAL_SPEC is',I
	END IF
!
	FIRST_TIME=.FALSE.
	FRAC_ELEC_HEATING=0.0D0
	FRAC_ION_HEATING=0.0D0
	FRAC_EXCITE_HEATING=0.0D0
	YE=0.0D0
	dXKT_MIN=(XKT_MAX-XKT_MIN)/(NKT-1)
	IF(METHOD .EQ. 'log')dXKT_MIN=XKT(2)-XKT(1)
!
! Set how the lectrons are ejected at high enegry. This will be the same for all
! depths. We can either inject at EMAX, or use a Bell-shaped curve near EMAX.
! E_INIT will be a normalisation constant to yield fractions for ionization etc.
!
	WRITE(LU_TH,*)'Constructing SOURCE'
	SOURCE=0.0D0
	IF (INJECT_DIRAC_AT_EMAX) THEN ! zero otherwise
	  E_INIT = XKT_MAX
	  WRITE(LU_TH,*)'Using Dirac delta function injection'
	ELSE
	  T1 = 0.0D0
	  T2 = DELTA_ENR_SOURCE
       	  IF (SOURCE_TYPE .EQ. 'CONSTANT')THEN
	    DO IKT=1,NKT
	      IF (XKT(IKT) .GE. (XKT_MAX-T2)) THEN
	        SOURCE(IKT) = 1.0D0
	        T1 = T1 + SOURCE(IKT)*dXKT(IKT)
              END IF
	    END DO
	    WRITE(LU_TH,*)'Using constant for electron injection'
	  ELSE IF(SOURCE_TYPE .EQ. 'BELL_SHAPE')THEN
	    DO IKT=1,NKT
	      IF (XKT(IKT) .GE. (XKT_MAX-T2)) THEN
                SOURCE(IKT) = 6.0D0 / T2**3 * (T2**2/4.D0 - (XKT(IKT)-XKT_MAX+T2/2.0D0)**2)
	        T1 = T1 + SOURCE(IKT)*dXKT(IKT)
              END IF
	    END DO
	    WRITE(LU_TH,*)'Using BELL_SHAPE for electron injection'
          END IF
	  SOURCE(:) = SOURCE(:) / T1
!
	  T1 = 0.0D0
	  T2 = 0.0D0
	  DO IKT=1,NKT
	    T1 = T1 + SOURCE(IKT)*dXKT(IKT)
	    T2 = T2 + XKT(IKT)*SOURCE(IKT)*dXKT(IKT)
	  END DO
	  E_INIT = T2
	END IF
!
! Write a summary of IONS, ionization potentials, and thrshold cross-sections.
!
	WRITE(LU_TH,'(X,A,T12,3X,A,2X,A,4X,A,2X,A,4X,A)')'Ion ID','IT','IST','SIG(IST)','SIG(IST+1)','Ion Pot.'
	DO IT=1,NUM_THD
	  IF(THD(IT)%PRES)THEN
            XION_POT = THD(IT)%ION_POT
	    IF(XION_POT .LT. XKT(1))THEN
	      IST=1
	      SIG1=THD(IT)%CROSS_SEC(1); SIG2=0.0D0
	    ELSE IF(XION_POT .GT. XKT(NKT))THEN
	      IST=NKT+1
	      SIG1=0.0D0; SIG2=0.0D0
	    ELSE
	      IST=GET_INDX_DP(XION_POT,XKT,NKT)
	      SIG1=THD(IT)%CROSS_SEC(IST); SIG2=THD(IT)%CROSS_SEC(IST+1)
	    END IF
	    ID=THD(IT)%LNK_TO_ION
	    WRITE(LU_TH,'(X,A,T12,2I5,2ES12.2,F12.2)')ION_ID(ID),IT,IST,SIG1,SIG2,XION_POT
	  END IF
	END DO
!
!
!	open(unit=lu_bethe,file='bethe_cross_chk',status='unknown')
	do_and_1st_time(:)=l_false
	DO DPTH_INDX=1,ND
!
	  WRITE(LU_TH,'(/,X,A,I3)')'Starting depth index: ',DPTH_INDX
!
! Decide on which species we will handle when computing the electron distribution.
! For each route, we sum up over all levels of that term. This approach will
! be independent of the SL assignments.
!
	  DO IT=1,NUM_THD
	    THD(IT)%N_ATOM=0.0D0
	    IF(THD(IT)%PRES)THEN
	      ID=THD(IT)%LNK_TO_ION
	      DO J=1,THD(IT)%N_STATES
	        K=THD(IT)%ATOM_STATES(J)
	        THD(IT)%N_ATOM=THD(IT)%N_ATOM+ATM(ID)%XzV_F(K,DPTH_INDX)
	      END DO
	    END IF
	    POP_VEC(IT)=THD(IT)%N_ATOM
	  END DO
!
! DO_THIS_ION and THD()%DO_THIS_ION_ROUTE contain similar information.
!   DO_THIS_ION is used for deciding by ION
!   DO_THIS_ION_ROUTE is used for deciding by the ionization route.
!
	  WRITE(LU_TH,'(X,A)')'Determining the top 10 ionization routes at this depth'
	  WRITE(LU_TH,'(2X,A,4X,A,7X,A,9X,A)')'IT','ID','Ion','Pop'
	  CALL INDEXX(NUM_THD,POP_VEC,INDX,L_FALSE)
	  DO_THIS_ION(:)=L_FALSE
	  THD(:)%DO_THIS_ION_ROUTE=L_FALSE
	  DO IT=1,NUM_THD
	    IF(POP_VEC(INDX(IT)) .LT. NT_OMIT_ION_SCALE*POP_VEC(INDX(1)) .AND. IT .GT. 10)EXIT
	    ID=THD(INDX(IT))%LNK_TO_ION
	    THD(INDX(IT))%DO_THIS_ION_ROUTE=.TRUE.
	    DO_THIS_ION(ID)=.TRUE.
	    do_and_1st_time(id)=.true.
	    WRITE(LU_TH,'(I4,2X,I4,4X,A6,ES12.3)')IT,ID,TRIM(ION_ID(ID)),POP_VEC(INDX(IT))
	  END DO
!
! Building the upper-diagonal matrix mat(1:nkt,1:nkt)
!
	  WRITE(LU_TH,*)'Constructing MAT'
!
! Get the energy losses through Coulomb interaction
!
	  WRITE(LU_TH,*)'Electron density=',ED(DPTH_INDX)  !ED(DPTH_INDX)
	  CALL GET_LELEC(LELEC,XKT,NKT,ED(DPTH_INDX))  !,ED(DPTH_INDX))
	  IF(VERBOSE_OUTPUT)THEN
	    DO IKT=1,NKT
	      WRITE(100,*)IKT,XKT(IKT),LELEC(IKT)
	    END DO
	  END IF
!
! Initialize MAT, and the add Diagonal term for Coulomb interaction with 
! thermal electrons.
!
	  MAT=0.0D0
	  DO IKT=1,NKT
	    MAT(IKT,IKT) = LELEC(IKT)
	  END DO
!
!
	  CALL TUNE(1,'MAT_ION')
	  IF(INCLUDE_IONIZATION)THEN
!
! Loop over all species and ionization stages
!
	    DO IT=1,NUM_THD
	      IF(THD(IT)%PRES .AND. THD(IT)%DO_THIS_ION_ROUTE)THEN
	        NATOM=THD(IT)%N_ATOM
                XION_POT = THD(IT)%ION_POT
	        WRITE(LU_TH,'(A,I3,A,F7.2,A)')' Ionization potential for IT=',IT,' is ',XION_POT,' eV'
!
	        IF(XION_POT .LT. XKT(1))THEN
	           IST=1
	        ELSE IF(XION_POT .GT. XKT(NKT))THEN
	           IST=NKT+1
	        ELSE
	           IST=GET_INDX_DP(XION_POT,XKT,NKT)
	        END IF
! 
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(IKT,IKTP,EKT,EKTP,EMIN,EMAX,XCROSS,NA_dX_CROSS)
	        DO IKTP=IST,NKT
                  XCROSS = THD(IT)%CROSS_SEC(IKTP)
	          IF(XCROSS .GT. 0.0D0)THEN
	            EKTP = XKT(IKTP)
                    EMAX = MIN(0.5d0*(EKTP+XION_POT),XKT(NKT))
	            NA_dX_CROSS=NATOM*dXKT(IKTP)*XCROSS
!
	            DO IKT=1,IKTP
	              EKT = XKT(IKT)
	              EMIN = MAX(EKTP-EKT,XION_POT)
	              IF(EMAX .GT. EMIN)THEN
	                MAT(IKT,IKTP) = MAT(IKT,IKTP) + NA_dX_CROSS*INTSIGC(EKTP,XION_POT,EMIN,EMAX)
	              END IF
!
                      IF (EKTP .GT. (2.0D0*EKT+XION_POT)) THEN
                        EMIN = EKT+XION_POT
                        MAT(IKT,IKTP) = MAT(IKT,IKTP) - NA_dX_CROSS*INTSIGC(EKTP,XION_POT,EMIN,EMAX)
                      END IF
	            END DO
	          END IF
!
	        END DO
!$OMP END PARALLEL DO
	      END IF    !Species present
	    END DO	!Looping over ionization stage / shell
	  END IF
	  CALL TUNE(2,'MAT_ION')
!
! Collisional excitations.
!
	CALL TUNE(1,'MAT_EXCITE')
	IF (INCLUDE_EXCITATION)THEN
	  WRITE(LU_TH,*)'Beginning excitation section'
!
! We only do excitations for ions with a fractional population > NT_OMIT_ION_SCALE.
! We exclude the final ion when computing the species and fractional population.
! We do at least one ionization stage for each species.
!
	  ION_SUM(:)=0.0D0
	  DO_THIS_ION_EXC(:)=.FALSE.
	  DO ISPEC=1,NUM_SPECIES
	    T1=0.0D0
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	     IF(ATM(ID)%XzV_PRES)ION_SUM(ID)=SUM(ATM(ID)%XzV_F(:,DPTH_INDX))
	     T1=T1+ION_SUM(ID)
	    END DO
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      IF(ION_SUM(ID)/T1 .GT. NT_OMIT_ION_SCALE)DO_THIS_ION_EXC(ID)=.TRUE.
	      SPEC_SUM(ID)=T1
	    END DO
	  END DO
!
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. DO_THIS_ION_EXC(ID))THEN
!
! Determine levels from which excitation will occur.
!
	      T1=ION_SUM(ID)
	      MAX_LOW_LEV=ATM(ID)%NXzV_F
	      DO I=ATM(ID)%NXzV_F,1,-1
	        IF(ATM(ID)%XzV_F(I,DPTH_INDX)/SPEC_SUM(ID) .GT. NT_OMIT_LEV_SCALE)EXIT
	        MAX_LOW_LEV=I
	      END DO
	      NLOW_MAXS(ID)=MAX(NLOW_MAXS(ID),MAX_LOW_LEV)
!
	      DO I=1,MAX_LOW_LEV
	        NL=I
	        J=I+1
	        DO WHILE(J .LE. ATM(ID)%NXzV_F)
	          NUP=J
!
! Now weight dE by Qnn rather than by AXzV_F(NL,NUP). Avoids division by zero when
! oscilator strength is zero.
!
	          CALL BETHE_APPROX_V5(Qnn,NL,NUP,XKT,dXKT_ON_XKT,NKT,ID,DPTH_INDX)
	          IF(Qnn(NKT) .GT. 0.0D0)THEN
	            dE=Qnn(NKT)*Hz_TO_eV*(ATM(ID)%EDGEXZV_F(NL)-ATM(ID)%EDGEXZV_F(J))
	            QSUM=Qnn(NKT)
	            K=J
	            DO WHILE(K+1 .LE. ATM(ID)%NXzV_F)
	              IF(Hz_TO_eV*(ATM(ID)%EDGEXZV_F(K+1)-ATM(ID)%EDGEXZV_F(J)) .GE. 1.0)EXIT
	              K=K+1
	              NUP=K
	              CALL BETHE_APPROX_V5(Qnn_TMP,NL,NUP,XKT,dXKT_ON_XKT,NKT,ID,DPTH_INDX)
	              IF(Qnn_TMP(NKT) .GT. 0.0D0)THEN
	                CALL PAR_VEC_SUM(QNN,Qnn_TMP,NKT) 
	                dE=dE+Qnn_TMP(NKT)*Hz_TO_eV*(ATM(ID)%EDGEXZV_F(NL)-ATM(ID)%EDGEXZV_F(K))
	                QSUM=QSUM+Qnn_TMP(NKT)
	              END IF
	            END DO
	            dE=dE/QSUM			!Weighted average energy of transition array.
	            J=K
!
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(IKT,IKT0,IKTP,IKTPN,EKT,EMAX)
	            DO IKT=1,NKT
	              EKT = XKT(IKT)
	              IKT0 = IKT
!
! Need to find the upper bound iktpn for which ektp = ekt + Excitation_Energy
!
	              EMAX = EKT+dE
	              IKTPN=IKT0 
	              IF (EMAX .GT. XKT(NKT-1)) then
	                IKTPN = NKT
	              ELSE
	                IF(XKT_METHOD .EQ. 'lin')THEN
	                  IKTPN=NINT( (EMAX-EKT)/dXKT_MIN+IKT0 )
	                ELSE
	                  DO WHILE (XKT(IKTPN) .LT. EMAX)
	                    IKTPN = IKTPN + 1
	                  END DO
	                END IF
	              END IF
!
	              DO IKTP=IKT0,IKTPN
	                MAT(IKT0,IKTP) = MAT(IKT0,IKTP) + Qnn(IKTP)
	              END DO
!
	            END DO	!Loop over energy
!$OMP END PARALLEL DO
!
	          END IF        !Check if f(i,j)=0
	          J=J+1
	        END DO		!Loop over upper level
	      END DO		!Loop over lower level
	    END IF		!Ionization stage present
	  END DO		!Loop over ionization stage
!	  close(lu_bethe)
	END IF			!Include excitations?
	CALL TUNE(2,'MAT_EXCITE')
!
! Check no diagonal term of mat is zero
!
	 DO IKT=1,NKT
	   IF (MAT(IKT,IKT) .EQ. 0.0D0) THEN
	     WRITE(LU_ER,*) 'Zero diagonal element in mat ',IKT,SOURCE(IKT)
	     STOP
	   END IF
	 END DO
!
! We now build the right-hand side and the matrix corresponding to
! Eq.7 of KF92.
!
	 IF (INJECT_DIRAC_AT_EMAX) THEN
	   RHS(1:NKT) = 1.0D0
	 ELSE
	   RHS(NKT) = 0.0D0
	   IF (SOURCE_TYPE .EQ. 'CONSTANT') RHS(NKT) = SOURCE(NKT) * dXKT(NKT)
	   DO IKT=NKT-1,1,-1
	     RHS(IKT) = RHS(IKT+1) + SOURCE(IKT) * dXKT(IKT)
	   END DO
	 END IF
!
! Solve for the degradation function YE
!
	 YE(NKT,DPTH_INDX)=0.0D0
	 IF(INJECT_DIRAC_AT_EMAX)YE(NKT,DPTH_INDX)=1.0D0/LELEC(NKT)
	 DO IKT=NKT-1,1,-1
	   T1 = 0.0D0
	   DO IKTP=IKT+1,NKT
	      T1 = T1 + MAT(IKT,IKTP)*YE(IKTP,DPTH_INDX)
	   END DO
	   YE(IKT,DPTH_INDX) = (RHS(IKT) - T1) / MAT(IKT,IKT) ! diagonal term is never zero (lelec contrib.)
	 END DO
!
! Compute electron heating
!
	    DO IKT=1,NKT
	      FRAC_ELEC_HEATING(DPTH_INDX)=FRAC_ELEC_HEATING(DPTH_INDX) + YE(IKT,DPTH_INDX)*LELEC(IKT)*dXKT(IKT)
	    END DO
	    T1=XKT(1)*LELEC(1)*YE(1,DPTH_INDX)/E_INIT
	    FRAC_ELEC_HEATING(DPTH_INDX)=FRAC_ELEC_HEATING(DPTH_INDX)/E_INIT
	    FRAC_ELEC_HEATING(DPTH_INDX)=FRAC_ELEC_HEATING(DPTH_INDX)+T1
	    WRITE(LU_TH,'(A,2ES12.4)')'(Fraction into electron heating and correction are',FRAC_ELEC_HEATING(DPTH_INDX),T1
!
!	  END DO		!Depth index
!
! Compute the ionization and heating contributions
!
	    IF(INCLUDE_IONIZATION .AND. DO_LOC_CHK)THEN
	      CALL TUNE(1,'CHK_ION')
	      WRITE(LU_TH,*)'Checking ionization energy'
	      WRITE(LU_TH,'(4X,A,2X,A,4X,A,2X,A,3X,A)')'IT','IST','SIG(IST)','SIG(IST+1)','dE(ION)/E'
	      DO IT=1,NUM_THD
	        IF(THD(IT)%PRES .AND. THD(IT)%DO_THIS_ION_ROUTE)THEN
	          NATOM=THD(IT)%N_ATOM
                  XION_POT = THD(IT)%ION_POT
!
	          IF(XION_POT .LT. XKT(1))THEN
	             IST=1
	             SIG1=THD(IT)%CROSS_SEC(1); SIG2=0.0D0
	          ELSE IF(XION_POT .GT. XKT(NKT))THEN
	             IST=NKT+1
	             SIG1=0.0D0; SIG2=0.0D0
	          ELSE
	             IST=GET_INDX_DP(XION_POT,XKT,NKT)
	             SIG1=THD(IT)%CROSS_SEC(IST); SIG2=THD(IT)%CROSS_SEC(IST+1)
	          END IF
	          T1 = 0.0D0
	          DO IKTP=IST,NKT
	            XCROSS = THD(IT)%CROSS_SEC(IKTP)
	            DETAI_DE(IKTP) = DETAI_DE(IKTP) + &
	                NATOM * XION_POT / E_INIT * YE(IKTP,DPTH_INDX)*XCROSS
	            T1 = T1 + YE(IKTP,DPTH_INDX)*dXKT(IKTP)*XCROSS
	          END DO
	          WRITE(LU_TH,'(X,2I5,3ES12.2)')IT,IST,SIG1,SIG2,NATOM * XION_POT * T1 / E_INIT
	          FRAC_ION_HEATING(DPTH_INDX) = FRAC_ION_HEATING(DPTH_INDX) + NATOM * XION_POT * T1 / E_INIT
	        END IF
	      END DO		!Ionization route
	    CALL TUNE(2,'CHK_ION')
	   END IF
	END DO		!Depth indx
!
	IF(INCLUDE_EXCITATION .AND. DO_LOC_CHK)THEN
	  CALL TUNE(1,'CHK_EXCITE')
!
!$OMP PARALLEL PRIVATE( RATE,dE,T1,SCALER_SPEC_SUM,SCALER_ION_SUM,ISPEC,ID,I,J,NL,NUP,IKT,DPTH_INDX,MAX_LOW_LEV )
!$OMP DO SCHEDULE(DYNAMIC)
!
	  DO DPTH_INDX=1,ND
	    DO ISPEC=1,NUM_SPECIES
	      SCALER_SPEC_SUM=0.0D0
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        T1=SUM(ATM(ID)%XzV_F(:,DPTH_INDX))
	        SCALER_SPEC_SUM=SCALER_SPEC_SUM+T1
	      END DO
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        SCALER_ION_SUM=SUM(ATM(ID)%XzV_F(:,DPTH_INDX))/SCALER_SPEC_SUM
	        IF(SCALER_ION_SUM .GT. NT_OMIT_ION_SCALE)THEN
!
	          T1=SCALER_ION_SUM
	          MAX_LOW_LEV=ATM(ID)%NXzV_F
	          DO I=ATM(ID)%NXzV_F,1,-1
	            IF(ATM(ID)%XzV_F(I,DPTH_INDX)/SCALER_SPEC_SUM .GT. NT_OMIT_LEV_SCALE)EXIT
	            MAX_LOW_LEV=I
	          END DO
!
	          DO I=1,MAX_LOW_LEV
	            DO J=I+1,ATM(ID)%NXzV_F
	              NL=I; NUP=J
	              CALL TOTAL_BETHE_RATE_V3(RATE,NL,NUP,YE,XKT,dXKT,NKT,ID,DPTH_INDX,ND)
	              dE=Hz_TO_eV*(ATM(ID)%EDGEXZV_F(NL)-ATM(ID)%EDGEXZV_F(NUP))
	              FRAC_EXCITE_HEATING(DPTH_INDX)=FRAC_EXCITE_HEATING(DPTH_INDX)+dE*RATE*ATM(ID)%XzV_F(NL,DPTH_INDX)
	            END DO
	          END DO
	        END IF		!Do this ion
	      END DO		!Loop over ionization stage
	    END DO              !Loop over species
	    FRAC_EXCITE_HEATING(DPTH_INDX)=FRAC_EXCITE_HEATING(DPTH_INDX)/E_INIT
	  END DO		!Loop over depth
!$OMP END DO
!$OMP END PARALLEL
	    CALL TUNE(2,'CHK_EXCITE')
	END IF		!Include excitation and DO_LOC_CHK?
!
! Provide a depth dependent summary
!
	IF(DO_LOC_CHK)THEN
	  WRITE(LU_TH,'(/,X,A,/)')'Summary of energy channels'
	  WRITE(LU_TH,'(4(5X,A))'),' dE(ELEC)/E','  dE(ION)/E','  dE(EXC)/E','dE(TOTAL)/E'
	  DO I=1,ND
	    WRITE(LU_TH,'(4ES16.6)')FRAC_ELEC_HEATING(I),FRAC_ION_HEATING(I),FRAC_EXCITE_HEATING(I), &
	    FRAC_ELEC_HEATING(I)+FRAC_ION_HEATING(I)+FRAC_EXCITE_HEATING(I)
	  END DO
	END IF
!
! Below we compute the number of non thermal electrons arising from the input of non-thermal energy.
! This is the number on electrons per 1eV of energy input.
!
	WRITE(LU_TH,*)' '
	WRITE(LU_TH,*)'Integral of normalized YE as a function of depth'
	WRITE(LU_TH,'(1X,A,5X,A,10X,A)')'Depth','Int(YE)','Ne'
	DO DPTH_INDX=1,ND
	  T1=0.0D0
	  DO IKT=1,NKT
	   T1=T1+YE(IKT,DPTH_INDX)/SQRT(XKT(IKT))*dXKT(IKT)
	  END DO
	  T1=T1*SQRT(9.109389D-28/2.0D0/1.602177D-12)/E_INIT
	  WRITE(LU_TH,'(I6,2ES12.3)')DPTH_INDX,T1,ED(DPTH_INDX)
	END DO
!
! Close diagnostic file, forcing all output.
!
	CLOSE(UNIT=LU_TH)
!
! Before leaving the routine we normalize YE by the E_INIT. We can later scale by the actual energy input.
!
	YE=YE/E_INIT
!
! Output the degradation spectrum
!
	open(unit=lu_bethe,file='NON_THERM_DEGRADATION_SPEC',status='unknown')
	write(lu_bethe,'(1x,A,I5)')'Number of energy grid: ',NKT
	write(lu_bethe,'(1x,A,I5)')'Number of depth: ',ND
	write(lu_bethe,*)''
	write(lu_bethe,'(A)')'Energy grids (eV)'
	write(lu_bethe,'(1x,1p8e16.7)')(XKT(IKT),IKT=1,NKT)
	write(lu_bethe,*)''
	write(lu_bethe,'(A)')'Source'
	write(lu_bethe,'(1x,1p8e16.7)')(SOURCE(IKT),IKT=1,NKT)
	write(lu_bethe,*)''
	write(lu_bethe,'(A)')'Lelec'
	write(lu_bethe,'(1x,1p8e16.7)')(LELEC(IKT),IKT=1,NKT)
	write(lu_bethe,*)''
	do dpth_indx=1,nd
	  write(lu_bethe,'(A,I5)')'Depth: ',dpth_indx
	  write(lu_bethe,'(1x,1p8e16.7)')(YE(IKT,DPTH_INDX),IKT=1,NKT)
	  write(lu_bethe,*)''
	end do
	close(unit=lu_bethe)
!
	WRITE(LU_ER,*)'Exiting ELECTRON_NON_THERM_SPEC'
!
	RETURN
	END
