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
	REAL*8 WRK(ND)
	REAL*8 LOG_R(ND)
	REAL*8 DELTA_T_SECS
	REAL*8 OLD_SN_AGE_DAYS
	REAL*8 DEN_SCL_FAC
	REAL*8 T1,T2
	INTEGER NX
	INTEGER NSP
	INTEGER NISO
	INTEGER I,L,K
	INTEGER IS,IP
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER(LEN=200) STRING
	LOGICAL, SAVE :: FIRST=.TRUE.
	LOGICAL DONE
!
	LUER=ERROR_LU()
	WRITE(LUER,*)'Entering RD_SN_DATA'
	OPEN(UNIT=LU,FILE='SN_HYDRO_DATA',STATUS='OLD')
!
! Get the number of data points in the HYDRO model, and the number 
! of species.
!
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
	 WRITE(LUER,*)'Obtained non-POP vectors in RD_SN_DATA'
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
!	    WRITE(LUER,'(I5,A,I5,2ES14.4)')L,TRIM(ISO_SPEC_HYDRO(L)),BARY_HYDRO(L),ISO_HYDRO(1,L),ISO_HYDRO(NX,L)
	    STRING=' '
	  END DO
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
	  WRITE(LUER,*)'RD_SN_DATA has read ED and T'
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
	WRITE(LUER,*)'Read SN populations in RD_SN_DATA'
!
	DO IS=1,NUM_ISOTOPES
	  ISO(IS)%OLD_POP_DECAY=ISO(IS)%OLD_POP
	END DO
	DO L=1,NISO
	  DONE=.FALSE.
	  DO K=1,NUM_ISOTOPES
	    IF(ISO_SPEC_HYDRO(L) .EQ. ISO(K)%SPECIES .AND. 
	1               BARY_HYDRO(L) .EQ. ISO(K)%BARYON_NUMBER)THEN
!	       CALL MON_INTERP(ISO(K)%OLD_POP,ND,IONE,LOG_R,ND,ISO_HYDRO(1:NX,L),NX,LOG_R_HYDRO,NX)
	       CALL LIN_INTERP(LOG_R,ISO(K)%OLD_POP,ND,LOG_R_HYDRO,ISO_HYDRO(1:NX,L),NX)
	       DONE=.TRUE.
	       EXIT
	    END IF
	  END DO
	  IF(.NOT. DONE)THEN
	    WRITE(LUER,*)'Error in RD_SN_DATA: isotope data not recognized'
	    WRITE(LUER,*)L,ISO_SPEC_HYDRO(L),BARY_HYDRO(L)
	    STOP
	  END IF
	END DO
	WRITE(LUER,*)'Read SN isotope populations in RD_SN_DATA'
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
	WRITE(6,*)'Normalized mass fractions in RD_SN_DATA'
	WRITE(6,*)'Maximum normalization factor was',MAXVAL(WRK)
	WRITE(6,*)'Minimum normalization factor was',MINVAL(WRK)
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
	      IF(ISO(IS)%ISPEC .EQ. PAR(IP)%ISPEC)THEN
	        WRK=WRK+ISO(IS)%OLD_POP
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
	  WRITE(LUER,*)'Normalized isotope populations in RD_SN_DATA'
	ELSE
	  DO IP=1,NUM_PARENTS
	    WRK(1:ND)=0.0D0
	    DO IS=1,NUM_ISOTOPES
	      IF(ISO(IS)%ISPEC .EQ. PAR(IP)%ISPEC)THEN
	        WRK=WRK+ISO(IS)%OLD_POP
	      END IF
	    END DO
	    T1=1.0D-100
	    DO I=1,ND
	       T2=(POP_SPECIES(I,PAR(IP)%ISPEC)+T1)/(WRK(I)+T1)-1.0D0
	       IF(SPECIES_PRES(PAR(IP)%ISPEC) .AND. ABS(T2) .GT. 0.20D0)THEN
	         WRITE(LUER,*)'Error in RD_SN_DATA: inconsistent total/isotope populations'
	         WRITE(LUER,*)'Species:',SPECIES(PAR(IP)%ISPEC)
	         WRITE(LUER,*)'Depth=',I,'  Fractional difference=',T2
	         STOP
	       END IF
	    END DO 
	    POP_SPECIES(:,PAR(IP)%ISPEC)=WRK
	  END DO
	END IF
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
! Correct populations for radioactive decays.
!
	IF(SN_AGE_DAYS.LT. OLD_SN_AGE_DAYS)THEN
	  WRITE(LUER,*)'Error in RD_SN_DATA'
	  WRITE(LUER,*)'Invalid epochs'
	  WRITE(LUER,*)SN_AGE_DAYS,OLD_SN_AGE_DAYS
	  STOP
	END IF
	DELTA_T_SECS=24.0D0*3600.0D0*(SN_AGE_DAYS-OLD_SN_AGE_DAYS)
	CALL DO_SPECIES_DECAYS(DELTA_T_SECS,ND)
        DO IS=1,NUM_ISOTOPES
          ISO(IS)%POP=ISO(IS)%OLD_POP_DECAY
        END DO
	DO IP=1,NUM_PARENTS
	  POP_SPECIES(1:ND,PAR(IP)%ISPEC)=PAR(IP)%OLD_POP_DECAY
	END DO
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
	CALL GET_NON_LOCAL_GAMMA_ENERGY(V,ND,LU)
!
	CALL OUT_SN_POPS_V3('SN_DATA_INPUT_CHK',SN_AGE_DAYS,USE_OLD_MF_OUTPUT,ND,LU)
!
	DEALLOCATE (R_HYDRO, LOG_R_HYDRO, V_HYDRO, SIGMA_HYDRO, T_HYDRO, DENSITY_HYDRO )
	DEALLOCATE (WRK_HYDRO, SPEC_HYDRO, POP_HYDRO, ATOM_DEN_HYDRO, ELEC_DEN_HYDRO)
	WRITE(LUER,*)'Exiting RD_SN_DATA'
!
	RETURN
	END
