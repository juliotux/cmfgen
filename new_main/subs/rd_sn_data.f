!
! Subroutine to read in hydrodynamical data data from a SN model.
! Routine reads in T, Density, and mass-fractions, and interpolates
! them onto the CMFGEN grid.
!
	SUBROUTINE RD_SN_DATA(ND,NEW_MODEL,LU)
	USE CONTROL_VARIABLE_MOD
	USE NUC_ISO_MOD
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered: 31-Aug-2016 : Better error reporting -- we check taht all chain isotopes are present.
! Altered: 15-May-2016 : Fixed minor bug when checking size of isotope abundance changes.
! Altered: 01-Mar-2016 : Changed to allow handling of a standard NUC_DECAY_DATA file.
!                         Code checks availability of decay route. This is important
!                         when a species but not isotopes are included [17-Feb-2016].
!                         Still may need to alter NUC_DECAY_DATA file is switch in middle of
!                         a model sequence (since a species with a short life may suddenly
!                         start to decay). A warning is ouput in such cases.
! Altered 15-Feb-2016: ISOTOPE_COUNTER added so that the same NUC_DECAY_DATA file can be used
!                         for all species.
! Altered 30-Jan-2015: Mass of each isotope is now output to SPECIES_MASSES
!                        Now call GET_NON_LOCAL_GAMMA_ENERGY_V2
! Altered  8-Nov-2009: Isotope/total population consistency check included.
! Altered 24-Jun-2009: PURE_HUBBLE_FLOW now stored in CONTROL_VARIABLE_MOD
! Altered 02-Apr-2009: Now set fixed abundance specifiers, normally read in, to surface values.
! Altered 11-Feb-2009: Use isotope data, when available, to atomic number fractions.
! Altered 28-Jan-2009: If inner & outer radii differ from passed R values (via MOD_CMFGEN)
!                        by less than 1 part in 10^8, they are made identical. This removes
!                        a problem with older models where a slightly different technique
!                        is used to compute the R grid.
! Altered  1-Jan-2009: Density scaling now valid for an arbitrary (but time independent
!                        in Lagrangian frame) velocity law.
! Altered 12-Feb-2008: Fixed pop in calculation of POP_ATOM (only counts species present).
!                      Now set species not present to have zero population.
!
	INTEGER LU
	INTEGER ND
	LOGICAL NEW_MODEL
!
! Local arrays & variables.
!
	REAL*8, ALLOCATABLE :: R_HYDRO(:)
	REAL*8, ALLOCATABLE :: LOG_R_HYDRO(:)
	REAL*8, ALLOCATABLE :: V_HYDRO(:)
	REAL*8, ALLOCATABLE :: SIGMA_HYDRO(:)
	REAL*8, ALLOCATABLE :: T_HYDRO(:)
	REAL*8, ALLOCATABLE :: DENSITY_HYDRO(:)
	REAL*8, ALLOCATABLE :: ATOM_DEN_HYDRO(:)
	REAL*8, ALLOCATABLE :: ELEC_DEN_HYDRO(:)
	REAL*8, ALLOCATABLE :: POP_HYDRO(:,:)
	REAL*8, ALLOCATABLE :: ISO_HYDRO(:,:)
	REAL*8, ALLOCATABLE :: WRK_HYDRO(:)
	INTEGER, ALLOCATABLE :: BARY_HYDRO(:)
	CHARACTER(LEN=10), ALLOCATABLE :: SPEC_HYDRO(:)
	CHARACTER(LEN=10), ALLOCATABLE :: ISO_SPEC_HYDRO(:)
!
	REAL*8 MASS_SPECIES(100)
	REAL*8 WRK(ND)
	REAL*8 LOG_R(ND)
	REAL*8 DELTA_T_SECS
	REAL*8 OLD_SN_AGE_DAYS
	REAL*8 DEN_SCL_FAC
	REAL*8 T1,T2
	REAL*8 MF_IB,MF_OB
	REAL*8 ISO_MF_IB,ISO_MF_OB
!
	INTEGER NX
	INTEGER NSP
	INTEGER NISO
	INTEGER I,L,K,CNT
	INTEGER IS,IP,ID,IN
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER(LEN=200) STRING
	LOGICAL, SAVE :: FIRST=.TRUE.
	LOGICAL DONE
	LOGICAL FIRST_WARN
!
	LUER=ERROR_LU()
	WRITE(LUER,'(/,A)')' Entering RD_SN_DATA'
	OPEN(UNIT=LU,FILE='SN_HYDRO_DATA',STATUS='OLD',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening SN_HYDRO_DATA ins rd_sn_data.f'
	  WRITE(LUER,*)'IOSTAT=',IOS
	  STOP
	END IF
!
! Get the number of data points in the HYDRO model, and the number 
! of species.
!
	  NX=0; NSP=0; NISO=0; OLD_SN_AGE_DAYS=0.0D0
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
	    ELSE IF(INDEX(STRING,'Number of isotopes:') .NE. 0)THEN
	      I=INDEX(STRING,'isotopes:')+9
	      READ(STRING(I:),*)NISO
	    ELSE IF(INDEX(STRING,'Time(days) since explosion:') .NE. 0)THEN
	      I=INDEX(STRING,'explosion:')+10
	      READ(STRING(I:),*)OLD_SN_AGE_DAYS
	    ELSE IF(INDEX(STRING,'Radius grid') .NE. 0)THEN
	      IF(NSP .EQ. 0 .OR. NX .EQ. 0 .OR. NISO .EQ. 0 .OR. OLD_SN_AGE_DAYS.EQ. 0.0D0)THEN
	        WRITE(LUER,*)'Error in RD_SN_DATA'
	        WRITE(LUER,*)'ND, NSP, NISO or model time undefined'
	        STOP
	      ELSE
	        EXIT
	      END IF
   	    END IF
	  END DO
!
	  ALLOCATE (R_HYDRO(NX));           R_HYDRO=0.0D0
	  ALLOCATE (LOG_R_HYDRO(NX));       LOG_R_HYDRO=0.0D0
	  ALLOCATE (V_HYDRO(NX));           V_HYDRO=0.0D0
	  ALLOCATE (SIGMA_HYDRO(NX));       SIGMA_HYDRO=0.0D0
	  ALLOCATE (T_HYDRO(NX));           T_HYDRO=0.0D0
	  ALLOCATE (ELEC_DEN_HYDRO(NX));    ELEC_DEN_HYDRO=0.0D0
	  ALLOCATE (ATOM_DEN_HYDRO(NX));    ATOM_DEN_HYDRO=0.0D0
	  ALLOCATE (DENSITY_HYDRO(NX));     DENSITY_HYDRO=0.0D0
	  ALLOCATE (WRK_HYDRO(NX));         WRK_HYDRO=0.0D0
	  ALLOCATE (SPEC_HYDRO(NSP));       SPEC_HYDRO=' '
	  ALLOCATE (POP_HYDRO(NX,NSP));     POP_HYDRO=0.0D0
!
	  ALLOCATE (ISO_HYDRO(NX,NISO));    ISO_HYDRO=0.0D0
	  ALLOCATE (BARY_HYDRO(NISO));      BARY_HYDRO=0
	  ALLOCATE (ISO_SPEC_HYDRO(NISO));  ISO_SPEC_HYDRO=' '
!
! Get basic HYDRO grid vectors.
!
	 DO WHILE (1 .EQ. 1)
	    DO WHILE (STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	      READ(LU,'(A)')STRING 
	    END DO
	    IF(INDEX(STRING,'Radius grid') .NE. 0)THEN
	      READ(LU,*)R_HYDRO
	    ELSE IF(INDEX(STRING,'Velocity') .NE. 0)THEN
	      READ(LU,*)V_HYDRO
	    ELSE IF(INDEX(STRING,'Sigma') .NE. 0)THEN
	      READ(LU,*)SIGMA_HYDRO
	    ELSE IF(INDEX(STRING,'Temperature') .NE. 0)THEN
	      READ(LU,*)T_HYDRO
	    ELSE IF(INDEX(STRING,'Density') .NE. 0)THEN
	      READ(LU,*)DENSITY_HYDRO
	    ELSE IF(INDEX(STRING,'Electron density') .NE. 0)THEN
	      READ(LU,*)ELEC_DEN_HYDRO
	    ELSE IF(INDEX(STRING,'Atom density') .NE. 0)THEN
	      READ(LU,*)ATOM_DEN_HYDRO
	    ELSE IF(INDEX(STRING,'ass fraction') .NE. 0)THEN
	      IF(R_HYDRO(1) .EQ. 0.0D0 .OR. 
	1	      T_HYDRO(1) .EQ. 0.0D0 .OR.
	1	      ATOM_DEN_HYDRO(1) .EQ. 0.0D0 .OR.
	1	      ELEC_DEN_HYDRO(1) .EQ. 0.0D0)THEN
	        WRITE(LUER,*)'Error reading SN data'
	        WRITE(LUER,*)'R, or T is zero'
	        STOP
	       ELSE
	         EXIT
	       END IF
	   END IF
	   STRING=' '
	 END DO
	 WRITE(LUER,*)'   Obtained non-POP vectors in RD_SN_DATA'
!
! We can now read in the mass-fractions.
!
	  DO L=1,NSP
	    DO WHILE( INDEX(STRING,'mass fraction') .EQ. 0)
	      READ(LU,'(A)')STRING 
	    END DO
	    STRING=ADJUSTL(STRING)
	    SPEC_HYDRO(L)=STRING(1:INDEX(STRING,' '))
	    READ(LU,*)(POP_HYDRO(I,L), I=1,NX)
	    STRING=' '
	  END DO
!
! We can now read in the isotope mass-fractions.
!
	  DO L=1,NISO
	    DO WHILE( INDEX(STRING,'mass fraction') .EQ. 0)
	      READ(LU,'(A)')STRING 
	    END DO
	    STRING=ADJUSTL(STRING)
	    I=INDEX(STRING,' ')
	    ISO_SPEC_HYDRO(L)=STRING(1:I-1)
	    READ(STRING(I:),*)BARY_HYDRO(L)
	    READ(LU,*)(ISO_HYDRO(I,L), I=1,NX)
	    STRING=' '
	  END DO
!
! Check no more isotopes
!
	  STRING=' '
	  DO WHILE( INDEX(STRING,'mass fraction') .EQ. 0)
	     READ(LU,'(A)',END=1000)STRING 
	  END DO
	  WRITE(LUER,*)'Error in RD_SN_DATA -- more isotopic mass fractions in SN_HYDRO_DATA' 
	  WRITE(LUER,*)TRIM(STRING)
	  STOP
1000	  CONTINUE
	CLOSE(UNIT=LU)
!
	IF(PURE_HUBBLE_FLOW)THEN
	  T1=24.0D0*3600.0D0*1.0D+05*OLD_SN_AGE_DAYS/1.0D+10
	  DO I=1,NX
	    T2=V_HYDRO(I)
	    V_HYDRO(I)=R_HYDRO(I)/T1
	    WRITE(233,*)I,T2,V_HYDRO(I)
	    SIGMA_HYDRO(I)=0.0D0
	  END DO
	END IF
!
! Correct populations for SN expansion. We do not correct
! POP_HYDRO and ISO_HYDRO as these are mass-fractions.
!
	T1=24.0D0*3600.0D0*1.0D+05*(SN_AGE_DAYS-OLD_SN_AGE_DAYS)/1.0D+10
	DO I=1,NX
	  T2=(1.0D0+T1*V_HYDRO(I)/R_HYDRO(I)*(SIGMA_HYDRO(I)+1.0D0))*
	1                   (1.0D0+T1*V_HYDRO(I)/R_HYDRO(I))**2
	  R_HYDRO(I)=R_HYDRO(I)+T1*V_HYDRO(I)
	  DENSITY_HYDRO(I)=DENSITY_HYDRO(I)/T2
	  ATOM_DEN_HYDRO(I)=ATOM_DEN_HYDRO(I)/T2
	  ELEC_DEN_HYDRO(I)=ELEC_DEN_HYDRO(I)/T2
	END DO
	IF( ABS(R_HYDRO(1)/R(1)-1.0D0) .LE. 1.0D-07 )R_HYDRO(1)=R(1)
	IF( ABS(R_HYDRO(NX)/R(ND)-1.0D0) .LE. 1.0D-07 )R_HYDRO(NX)=R(ND)
	DO I=1,ND
	  VOL_EXP_FAC(I)=1.0D0/(1.0D0-T1*V(I)/R(I)*(SIGMA(I)+1.0D0))/
	1                   (1.0D0-T1*V(I)/R(I))**2
	END DO
!
! We can now interpolate from the HYDRO grid to the CMFGEN grid.
! Interpolations are don in R, in the log-log plane.
!
! If the model is not a NEW_MODEL, ED and T are set elsewhere, and thus we 
! don't want to corrupt them.
!
	LOG_R_HYDRO=LOG(R_HYDRO)
	LOG_R=LOG(R)
	IF(FIRST .AND. NEW_MODEL)THEN
	  T_HYDRO=LOG(T_HYDRO) 
	  CALL MON_INTERP(T,ND,IONE,LOG_R,ND,T_HYDRO,NX,LOG_R_HYDRO,NX)
	  T=EXP(T)
	  ELEC_DEN_HYDRO=LOG(ELEC_DEN_HYDRO) 
	  CALL MON_INTERP(ED,ND,IONE,LOG_R,ND,ELEC_DEN_HYDRO,NX,LOG_R_HYDRO,NX)
	  ED=EXP(ED)
	  WRITE(LUER,*)'   RD_SN_DATA has read ED and T'
	END IF
!
	WRK_HYDRO=LOG(DENSITY_HYDRO) 
	CALL MON_INTERP(DENSITY,ND,IONE,LOG_R,ND,WRK_HYDRO,NX,LOG_R_HYDRO,NX)
	DENSITY=EXP(DENSITY)
!
! Changed to LIN_INTERP as some models have humongous grid changes across
! grid points which can cause -ve values due to round-off errors.
! 
	POP_SPECIES=0.0D0
	DO L=1,NSP
	  DO K=1,NUM_SPECIES
	    IF(SPEC_HYDRO(L) .EQ. SPECIES(K))THEN
!	      CALL MON_INTERP(POP_SPECIES(1,K),ND,IONE,LOG_R,ND,POP_HYDRO(1,L),NX,LOG_R_HYDRO,NX)
	      CALL LIN_INTERP(LOG_R,POP_SPECIES(1,K),ND,LOG_R_HYDRO,POP_HYDRO(1,L),NX)
	   END IF
	  END DO
	END DO
	WRITE(LUER,*)'   Read SN populations in RD_SN_DATA'
!
	ISO(:)%READ_ISO_POPS=.FALSE.
	DO L=1,NISO
	  DONE=.FALSE.
	  DO K=1,NUM_ISOTOPES
	    IF(ISO_SPEC_HYDRO(L) .EQ. ISO(K)%SPECIES .AND. 
	1               BARY_HYDRO(L) .EQ. ISO(K)%BARYON_NUMBER)THEN
!	       CALL MON_INTERP(ISO(K)%OLD_POP,ND,IONE,LOG_R,ND,ISO_HYDRO(1:NX,L),NX,LOG_R_HYDRO,NX)
	       CALL LIN_INTERP(LOG_R,ISO(K)%OLD_POP,ND,LOG_R_HYDRO,ISO_HYDRO(1:NX,L),NX)
	       DONE=.TRUE.
	       ISO(K)%READ_ISO_POPS=.TRUE.
	       EXIT
	    END IF
	  END DO
	  IF(.NOT. DONE)THEN
	    WRITE(LUER,*)'Error in RD_SN_DATA: isotope data not recognized'
	    WRITE(LUER,*)L,ISO_SPEC_HYDRO(L),BARY_HYDRO(L)
	    STOP
	  END IF
	END DO
	WRITE(LUER,*)'   Read SN isotope populations in RD_SN_DATA'
!
	DO IS=1,NUM_ISOTOPES
	  IF(ISO(IS)%READ_ISO_POPS)THEN
	    DO I=1,ND
	     ISO(IS)%OLD_POP(I)=MAX(ISO(IS)%OLD_POP(I),MINIMUM_ISO_POP)
	    END DO
	  END IF
	END DO
!
! Ensure mass-fractions sum to unity. Two options: scale either
! the mass-fractions or the density. Here we scale the mass-fractions.
! To avoid possible confusion, we only allow for species explicitly
! included in the model.
!
	WRK(:)=0.0D0
	DO L=1,NUM_SPECIES
	  IF(SPECIES_PRES(L))THEN
	    WRK(:)=WRK(:)+POP_SPECIES(:,L)
	  ELSE
	    POP_SPECIES(:,L)=0.0D0
	  END IF
	END DO
	DO L=1,NUM_SPECIES
	  POP_SPECIES(:,L)=POP_SPECIES(:,L)/WRK(:)
	END DO
	WRITE(LUER,*)'   Normalized mass fractions in RD_SN_DATA'
	WRITE(LUER,*)'   Maximum normalization factor was',MAXVAL(WRK)
	WRITE(LUER,*)'   Minimum normalization factor was',MINVAL(WRK)
!
! Compute the mass of each species present in the ejecta.
!
	MASS_SPECIES=0.0D0
	DO L=1,NUM_SPECIES
	  IF(SPECIES_PRES(L))THEN
	    DO K=1,ND
	      WRK(K)=POP_SPECIES(K,L)*DENSITY(K)*R(K)*R(K)
	    END DO
	    CALL LUM_FROM_ETA(WRK,R,ND)
	    MASS_SPECIES(L)=4.0D0*3.1416D0*SUM(WRK(1:ND))/1.989D+03
	  END IF
	END DO
!
! Now compute the atomic population of each species.
!
	DO L=1,NUM_SPECIES
	  POP_SPECIES(:,L)=POP_SPECIES(:,L)*DENSITY(:)/AT_MASS(L)/1.66D-24
	END DO
	DO IS=1,NUM_ISOTOPES
	  ISO(IS)%OLD_POP=ISO(IS)%OLD_POP*DENSITY/ISO(IS)%MASS/1.66D-24
	END DO
!
! When isotopes are present, we need to use the individual mass-fractions to
! compute the correct number of species atoms that are present. In such cases
! the POP_PSECIES computed here overrides that computed above.
!
	IF(USE_OLD_MF_SCALING)THEN
!
! Ensure isotope populations sum exactly to total species population
! Isotope populations of species not-present get set to zero.
! Only necessary if using old scaling approach.
!
	  DO IP=1,NUM_PARENTS
	    WRK(1:ND)=0.0D0
	    DO IS=1,NUM_ISOTOPES
	      IF(ISO(IS)%ISPEC .EQ. PAR(IP)%ISPEC .AND. ISO(IS)%READ_ISO_POPS)THEN
	        WRK=WRK+ISO(IS)%OLD_POP
	        PAR(IP)%DECAY_CHAIN_AVAILABLE=.TRUE.
	      END IF
	    END DO
	    DO IS=1,NUM_ISOTOPES
	      IF(ISO(IS)%ISPEC .EQ. PAR(IP)%ISPEC)THEN
	        IF(SPECIES_PRES(ISO(IS)%ISPEC))THEN
	          ISO(IS)%OLD_POP=ISO(IS)%OLD_POP*(POP_SPECIES(:,ISO(IS)%ISPEC)/WRK)
	        ELSE
	          ISO(IS)%OLD_POP=0.0D0
	        END IF
	      END IF
	    END DO
	  END DO
	  WRITE(LUER,*)'   Normalized isotope populations in RD_SN_DATA'
	ELSE
	  DO IP=1,NUM_PARENTS
	    WRK(1:ND)=0.0D0
	    PAR(IP)%DECAY_CHAIN_AVAILABLE=.FALSE.
	    DO IS=1,NUM_ISOTOPES
	      IF(ISO(IS)%ISPEC .EQ. PAR(IP)%ISPEC .AND. ISO(IS)%READ_ISO_POPS)THEN
	        WRK=WRK+ISO(IS)%OLD_POP
	        PAR(IP)%DECAY_CHAIN_AVAILABLE=.TRUE.
	      END IF
	    END DO
	    IF(PAR(IP)%DECAY_CHAIN_AVAILABLE)THEN
	      T1=1.0D-100
	      DO I=1,ND
	         T2=(POP_SPECIES(I,PAR(IP)%ISPEC)+T1)/(WRK(I)+T1)-1.0D0
	         IF(SPECIES_PRES(PAR(IP)%ISPEC) .AND. ABS(T2) .GT. 0.20D0)THEN
	           WRITE(LUER,*)'Error in RD_SN_DATA: inconsistent total/isotope populations'
	           WRITE(LUER,*)'Species:',SPECIES(PAR(IP)%ISPEC)
	           WRITE(LUER,*)'Depth=',I,'  Fractional difference=',T2
	           WRITE(LUER,*)WRK(I),POP_SPECIES(I,PAR(IP)%ISPEC)
	           STOP
	         END IF
	      END DO 
	      POP_SPECIES(:,PAR(IP)%ISPEC)=WRK
	    END IF
	  END DO
	END IF
!
	FIRST_WARN=.TRUE.
	DO IP=1,NUM_PARENTS
	  IF(PAR(IP)%DECAY_CHAIN_AVAILABLE)THEN
	  ELSE
	    IF(FIRST_WARN)THEN
	      WRITE(LUER,*)' '
	      WRITE(LUER,*)'Warning: Possible error in reading SN_HYDRO_DATA with RD_SN_DATA.'
	      WRITE(LUER,*)'The following species have ISOTOPE data present in NUC_DECAY_DATA',
	1                        '     but have no isotope data in SN_HYDRO_DATA.'
	      WRITE(LUER,*)'You many need to add the appropriate isotopic data'
	      WRITE(LUER,*)SPECIES(PAR(IP)%ISPEC)
	      FIRST_WARN=.FALSE.
	    ELSE
	      WRITE(LUER,*)SPECIES(PAR(IP)%ISPEC)
	    END IF
	  END IF
	END DO
	WRITE(LUER,*)' '
!
! If population is zero at some depths, but species is present, we will
! set to small value.
!
	DO L=1,NUM_SPECIES
	  IF(SPECIES_PRES(L))THEN
	    DO K=1,ND
	      POP_SPECIES(K,L)=MAX(1.0D-20,POP_SPECIES(K,L))
	    END DO
	  END IF
	END DO
!
! Set for species with decay chain data, but which may not have isotopic
! data present.
!
	DO IP=1,NUM_PARENTS
	  PAR(IP)%OLD_POP=POP_SPECIES(:,PAR(IP)%ISPEC)
	  PAR(IP)%OLD_POP_DECAY= PAR(IP)%OLD_POP
	END DO
!
! Correct populations for radioactive decays.
!
	IF(SN_AGE_DAYS.LT. OLD_SN_AGE_DAYS)THEN
	  WRITE(LUER,*)'Error in RD_SN_DATA'
	  WRITE(LUER,*)'Invalid epochs'
	  WRITE(LUER,*)SN_AGE_DAYS,OLD_SN_AGE_DAYS
	  STOP
	END IF
	DELTA_T_SECS=24.0D0*3600.0D0*(SN_AGE_DAYS-OLD_SN_AGE_DAYS)
	CALL DO_SPECIES_DECAYS_V2(INSTANTANEOUS_ENERGY_DEPOSITION,DELTA_T_SECS,ND)
	IF(ABS(DELTA_T_SECS/SN_AGE_DAYS) .LT. 1.0D-05 .AND. 
	1     (INCL_DJDT_TERMS .OR. DO_CO_MOV_DDT) )THEN
	  WRITE(LUER,'(/,/,1X,A80)')('*',I=1,79)
	  WRITE(LUER,*)'Error in RD_SN_DATA'
	  WRITE(LUER,*)'You have requested time depndent calculations but have a very small time step'
	  WRITE(LUER,*)'DELTA_T_SECS=',DELTA_T_SECS
	  WRITE(LUER,*)'SN_AGE_DAYS=',SN_AGE_DAYS
	  WRITE(LUER,*)'OLD_SN_AGE_DAYS=',OLD_SN_AGE_DAYS
	  STOP
	END IF
!
! This checks that ISOTOPES for each ACTIVE chain is avaliable. An actve chain
! is one that has at least one isotope present.
!
	DO IN=1,NUM_DECAY_PATHS
	  IS=NUC(IN)%LNK_TO_ISO
	  ID=NUC(IN)%DAUGHTER_LNK_TO_ISO
	  IF(ISO(IS)%READ_ISO_POPS .NEQV. ISO(ID)%READ_ISO_POPS)THEN
	    WRITE(LUER,*)' '
	    WRITE(LUER,*)' Error -- you have read in one isotope in a decay chain but not both'
	    WRITE(LUER,'(1X,A5,I4,3X,A,:L1)') NUC(IN)%SPECIES,NUC(IN)%BARYON_NUMBER,'Read=',ISO(IS)%READ_ISO_POPS
	    WRITE(LUER,'(1X,A5,I4,3X,A,:L1)') NUC(IN)%DAUGHTER,NUC(IN)%DAUGHTER_BARYON_NUMBER,'Read=',ISO(ID)%READ_ISO_POPS
	    STOP
	  END IF
	END DO
!
        DO IS=1,NUM_ISOTOPES
          ISO(IS)%POP=ISO(IS)%OLD_POP_DECAY
        END DO
	DO IP=1,NUM_PARENTS
	  IF(PAR(IP)%DECAY_CHAIN_AVAILABLE)THEN
	    POP_SPECIES(1:ND,PAR(IP)%ISPEC)=PAR(IP)%OLD_POP_DECAY
	  END IF
	END DO
!
	FIRST_WARN=.TRUE.
	DO IS=1,NUM_ISOTOPES
	  IF(ISO(IS)%READ_ISO_POPS)THEN
	    T1=SUM(ISO(IS)%OLD_POP)
	    IF(T1 .NE. 0.0D0)THEN
	      T1=0.0D0
	      DO I=1,ND
	        T2=1.0D0-ISO(IS)%OLD_POP_DECAY(I)/ISO(IS)%OLD_POP(I)
	        IF( (T2 .LE. -2.0D0 .OR. T2 .GT. 0.67D0) .AND. ISO(IS)%OLD_POP(I) .GT. 1.0D-10)THEN
	          IF(FIRST_WARN)THEN
	            FIRST_WARN=.FALSE.
	            WRITE(LUER,*)' '
	            WRITE(LUER,*)'WARNING from RD_SN_DATA '
	            WRITE(LUER,*)'The following isotopes have changed their abundance by over a factor of 3.'
	            WRITE(LUER,*)'This may have occurred because you are using a revised NUCLEAR decay data file.'
	            WRITE(LUER,*)'Do diff SN_HYDRO_DATA SN_HYDRO_FOR_NEXT_MODEL to see changes.'
	            WRITE(LUER,*)'This may affect the heating, partyicularly when using energy'//
	1                        ' deposition averaged over time.'
	            WRITE(LUER,'(3X,A,T12,I3)')TRIM(ISO(IS)%SPECIES),ISO(IS)%BARYON_NUMBER
	          ELSE
	            WRITE(LUER,'(3X,A,T12,I3)')TRIM(ISO(IS)%SPECIES),ISO(IS)%BARYON_NUMBER
	          END IF
	          EXIT
	        END IF
	      END DO
	    END IF
	  END IF
	END DO
	IF(.NOT. FIRST_WARN)WRITE(LUER,*)' '
!
! Compute the mass for each isotope. We also output the mass fractions at the boundaries.
!
	OPEN(UNIT=LU,FILE='SPECIES_MASSES',STATUS='UNKNOWN')
	WRITE(LU,'(A,3X,A,5X,A,13X,8X,A,8X,A)')'Species ','M(Msun)','BN','MF(OB)','MF(IB)'
	DO L=1,NUM_SPECIES
	  IF(SPECIES_PRES(L))THEN
	    T2=0.0D0
	    CNT=0
!
! We use these MF's if no isotopes.
!
	    MF_OB=1.66D-24*AT_MASS(L)*POP_SPECIES(1,L)/DENSITY(1)
	    MF_IB=1.66D-24*AT_MASS(L)*POP_SPECIES(ND,L)/DENSITY(ND)
	    DO IS=1,NUM_ISOTOPES
	      IF(ISO(IS)%SPECIES .EQ. SPECIES(L))THEN
	        CNT=CNT+1
	        IF(CNT .EQ. 1)THEN
	          WRITE(LU,'(A)')' '
	          MF_OB=0.0D0; MF_IB=0.0D0
	        END IF
	        DO K=1,ND
	          WRK(K)=ISO(IS)%POP(K)*R(K)*R(K)
	        END DO
	        CALL LUM_FROM_ETA(WRK,R,ND)
	        T1=4.0D0*3.1416D0*1.66D-24*SUM(WRK(1:ND))*ISO(IS)%MASS/1.989D+03
	        ISO_MF_OB=ISO(IS)%POP(1)*1.66D-24*ISO(IS)%MASS/DENSITY(1)
	        ISO_MF_IB=ISO(IS)%POP(ND)*1.66D-24*ISO(IS)%MASS/DENSITY(ND)
	        MF_OB=MF_OB+ISO_MF_OB
	        MF_IB=MF_IB+ISO_MF_IB
	        IF(T1 .NE. 0 .AND. ISO(IS)%STABLE)THEN
	          WRITE(LU,'(A,T8,ES11.3,4X,I3,5X,A,2ES14.3)')TRIM(SPECIES(L)),
	1            T1,ISO(IS)%BARYON_NUMBER,'  Stable',ISO_MF_OB,ISO_MF_IB
	          T2=T2+T1
	        ELSE IF(T1 .NE. 0)THEN
	          WRITE(LU,'(A,T8,ES11.3,4X,I3,5X,A,2ES14.3)')TRIM(SPECIES(L)),
	1            T1,ISO(IS)%BARYON_NUMBER,'Unstable',ISO_MF_OB,ISO_MF_IB
	          T2=T2+T1
	        END IF
	      END IF
	    END DO
	    IF(T2 .NE. 0.0D0)MASS_SPECIES(L)=T2
	    WRITE(LU,'(A,T8,ES11.3,20X,2ES14.3)')TRIM(SPECIES(L)),MASS_SPECIES(L),
	1             MF_OB,MF_IB
	  END IF
	END DO
	WRITE(LU,'(/,A,ES12.4)')'Total ejecta mass of model is ',SUM(MASS_SPECIES)
	CLOSE(UNIT=LU)
!
! Compute total ATOM population.
!
	POP_ATOM(:)=0.0D0
	DO L=1,NUM_SPECIES
	  IF(SPECIES_PRES(L))THEN
	    POP_ATOM(:)=POP_ATOM(:)+POP_SPECIES(:,L)
	  END IF
	END DO
!
! Set fixed abundance specifiers, normally read in, to surface values.
!
	DO L=1,NUM_SPECIES
	  IF(SPECIES_PRES(L))AT_ABUND(L)=POP_SPECIES(1,L)/POP_ATOM(1)
	END DO
!
! Get non-local energy deposition, if important.
! This replaces that computed by DO_SPECIES_DECAYS computed earlier.
!
	CALL GET_NON_LOCAL_GAMMA_ENERGY_V2(R,V,ND,LU)
!
	CALL OUT_SN_POPS_V3('SN_DATA_INPUT_CHK',SN_AGE_DAYS,USE_OLD_MF_OUTPUT,ND,LU)
!
	DEALLOCATE (R_HYDRO, LOG_R_HYDRO, V_HYDRO, SIGMA_HYDRO, T_HYDRO, DENSITY_HYDRO )
	DEALLOCATE (ATOM_DEN_HYDRO, ELEC_DEN_HYDRO, POP_HYDRO, ISO_HYDRO, WRK_HYDRO)
	DEALLOCATE (BARY_HYDRO, SPEC_HYDRO, ISO_SPEC_HYDRO)
	WRITE(LUER,'(A,/)')' Exiting RD_SN_DATA'
!
	RETURN
	END
