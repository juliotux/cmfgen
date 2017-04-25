!
! Program to check energy conservation as a function of time for a CMFGEN SN model sequence.
! The spectra directory must be specified in the file DIRECTORIES.
!
! Program reads:
!              VADAT
!              RVTJ
!              OBSFLUX
!              JH_AT_CURRENT_TIME
! from the specified directory.
!
	PROGRAM CHECK_ENERGY_CONS
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 28-AUg-2016: Output obsereved luminosity
!                       Output better error diagnostics
!
	INTEGER, PARAMETER :: NMAX=200			!Maximum number of models in time sequence
	INTEGER, PARAMETER :: ND_MAX=200		!Maximum number of depth points (in any model)
	INTEGER, PARAMETER :: NCF_MAX=500000		!Maximum nnumber of frequncies (in any model)
	INTEGER, PARAMETER :: LUIN=10
	INTEGER, PARAMETER :: LUOUT=40
!
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: R_RVTJ(:)
	REAL*8, ALLOCATABLE :: RJ(:,:)
	REAL*8, ALLOCATABLE :: HFLUX(:,:)
	REAL*8, ALLOCATABLE :: NU(:)
	REAL*8, ALLOCATABLE :: OBS_FREQ(:)
	REAL*8, ALLOCATABLE :: OBS_FLUX(:)
!
	REAL*8 SN_AGE(NMAX)
	REAL*8 VINF(NMAX)
	INTEGER ND
!
	REAL*8 CMF_LUM(ND_MAX)
	REAL*8 MECH_LUM(ND_MAX)
	REAL*8 DECAY_LUM(ND_MAX)
	REAL*8 DR4J_LUM(ND_MAX)
	REAL*8 INTERN_LUM(ND_MAX)
!
	REAL*8 E_RAD(NMAX)
	REAL*8 L_OBS(NMAX)
	REAL*8 L_CMF(NMAX)
	REAL*8 L_MECH(NMAX)
	REAL*8 L_DECAY(NMAX)
	REAL*8 L_DR4J(NMAX)
	REAL*8 L_INTERN(NMAX)
	REAL*8 E_CONS(NMAX)
!
	REAL*8 LAMC
	REAL*8 T1,T2,T3,T4,T5
	REAL*8 PI
	REAL*8 LUM_CONV_FACTOR
	REAL*8 C_Mms
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
! For plotting.
!
	REAL*8 XV(NCF_MAX)
	REAL*8 YV(NCF_MAX)
!
	INTEGER NCF
	INTEGER NCF_CONT
	INTEGER NMOD
	INTEGER I,J,K,ML
	INTEGER REC_LENGTH
	INTEGER ST_REC
	INTEGER IOS
	INTEGER IMIN,IMAX
!
	LOGICAL DONL_CMF_LUM
	LOGICAL DONL_MECH_LUM
	LOGICAL DONL_DECAY_LUM
	LOGICAL DONL_DR4J_LUM
	LOGICAL DONL_INTERN_LUM
!
	CHARACTER(LEN=140) DIR_NAME(NMAX)
	CHARACTER(LEN=140) FILE_NAME
	CHARACTER(LEN=80) FILE_DATE
	CHARACTER(LEN=200) STRING
!
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	SN_AGE=0.0D0
	VINF=0.0D0
	ND=0
	C_MMS=1.0D-08*SPEED_OF_LIGHT()
	PI=4.0D0*ATAN(1.0D0)
!
	FILE_NAME='DIRECTORIES'
	WRITE(6,*)' '
	CALL GEN_IN(FILE_NAME,'File with list of directories')
	OPEN(UNIT=LUIN,FILE=FILE_NAME,STATUS='OLD',ACTION='READ')
	DO I=1,NMAX
	  READ(LUIN,'(A)',END=100)DIR_NAME(I)
	  NMOD=I
	  K=LEN_TRIM(DIR_NAME(I))
	  IF(DIR_NAME(I)(K:K) .NE. '/')DIR_NAME(I)(K+1:K+1)='/'
	END DO
100	CLOSE(UNIT=10)
	WRITE(6,*)'Number of directory names read is:',NMOD
	WRITE(6,*)' '
!
	OPEN(UNIT=LUOUT,FILE='TIME_INTEGRATED_ENERGY_CHK',STATUS='UNKNOWN')
	WRITE(LUOUT,'(/,A)')' L are in units of L(sun)'
	WRITE(LUOUT,'(A)')' E(rad) is in ergs'
	WRITE(LUOUT,'(A)')' L_INTERN refers to De/Dt etc'
	WRITE(LUOUT,'(A)')' L_DR4J refers Dr^4/dT term when L_MECH is zero.'
	WRITE(LUOUT,'(A,/)')' L_DR4J refers Dr^3/dT term when L_MECH is non zero.'
!
	WRITE(LUOUT,'(4X,A,7(6X,A))')'AGE(d)','   L_OBS','   L_CMF','  L_MECH','  L_DR4J',
	1                         'L_INTERN',' L_DECAY','   E_RAD'
	DO I=1,NMOD
!
          FILE_NAME=TRIM(DIR_NAME(I))//'RVTJ'
          ND=0
	  OPEN(UNIT=LUIN,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	    IF(IOS .EQ. 0)THEN
              DO WHILE(ND .EQ. 0)
                READ(LUIN,'(A)')STRING
                IF(INDEX(STRING,'ND:') .NE. 0)THEN
                  K=INDEX(STRING,':')+1
                  READ(STRING(K:),*)ND
                END IF
              END DO
	      IF(ALLOCATED(R_RVTJ))DEALLOCATE(R_RVTJ)
	      ALLOCATE (R_RVTJ(ND))
	      STRING=' '
	      DO WHILE(INDEX(STRING,'Radius (10^10 cm)') .EQ. 0)
	        READ(10,'(A)')STRING
	      END DO
	      READ(10,*)R_RVTJ	
	      CLOSE(UNIT=10)
	    ELSE
	      WRITE(6,*)' '
	      WRITE(6,*)'Unable to open ',TRIM(FILE_NAME)
	      WRITE(6,*)'Assuming this, and subsequent models, are unavailable'
	      NMOD=I-1
	      EXIT
	      WRITE(6,*)' '
	    END IF
!
	  DONL_CMF_LUM=.FALSE.
	  DONL_MECH_LUM=.FALSE.
	  DONL_DECAY_LUM=.FALSE.
	  DONL_DR4J_LUM=.FALSE.
	  DONL_INTERN_LUM=.FALSE.
          LUM_CONV_FACTOR=16*ATAN(1.0D0)*1.0D+15*1.0D-23*((3.0856D+21)**2)/3.826D+33
!
	  FILE_NAME=TRIM(DIR_NAME(I))//'OBSFLUX'
	  OPEN(UNIT=LUIN,FILE=FILE_NAME,STATUS='OLD',ACTION='READ')
	  DO WHILE(1 .EQ. 1)
	    READ(LUIN,'(A)',END=200)STRING
	    IF(INDEX(STRING,'Continuum Frequencies') .NE. 0)THEN
	      K=INDEX(STRING,'(')+1
	      IF(K .EQ. 1)THEN
	        WRITE(6,*)'Assuming maximum number of frequenices is ',NCF_MAX
	        ALLOCATE(OBS_FREQ(NCF_MAX)); OBS_FREQ=0.0D0
	        READ(LUIN,*,ERR=50)(OBS_FREQ(J),J=1,NCF_MAX)
50	        NCF=NCF_MAX
	        DO WHILE(OBS_FREQ(NCF) .EQ. 0.0D0)
	          NCF=NCF-1
	        END DO
	        BACKSPACE(LUIN)		!as may have read in obs. intensity record. 
	      ELSE
	        READ(STRING(K:),*)NCF
	        ALLOCATE(OBS_FREQ(NCF))
	        READ(LUIN,*)OBS_FREQ(1:NCF)
	      END IF
	      ALLOCATE(OBS_FLUX(NCF))
	    END IF
	    IF(INDEX(STRING,'Observed intensity') .NE. 0 .AND. .NOT. DONL_CMF_LUM)THEN
	      READ(LUIN,*)OBS_FLUX(1:NCF)
              L_OBS(I)=0.0D0
              DO K=2,NCF
                 L_OBS(I)=L_OBS(I)+(OBS_FREQ(K-1)-OBS_FREQ(K))*(OBS_FLUX(K)+OBS_FLUX(K+1))
              END DO
              L_OBS(I)=LUM_CONV_FACTOR*L_OBS(I)*0.5D0
	      DEALLOCATE(OBS_FREQ,OBS_FLUX)
	    END IF
	    IF(INDEX(STRING,'Luminosity') .NE. 0 .AND. .NOT. DONL_CMF_LUM)THEN
	      READ(LUIN,*)CMF_LUM(1:ND)
	      L_CMF(I)=CMF_LUM(1)
	      DONL_CMF_LUM=.TRUE.
	    END IF
	    IF(INDEX(STRING,'Mechanical Luminosity') .NE. 0  .AND. .NOT. DONL_MECH_LUM)THEN
	      READ(LUIN,*)MECH_LUM(1:ND)
	      L_MECH(I)=SUM(MECH_LUM(1:ND))
	      DONL_MECH_LUM=.TRUE.
	    END IF
	    IF(INDEX(STRING,'Internal') .NE. 0 .AND. .NOT. DONL_INTERN_LUM)THEN
	      READ(LUIN,*)INTERN_LUM(1:ND)
	      L_INTERN(I)=SUM(INTERN_LUM(1:ND))
	      DONL_INTERN_LUM=.TRUE.
	    END IF
	    IF(INDEX(STRING,'Flux arrising') .NE. 0 .AND. .NOT. DONL_DR4J_LUM)THEN
	      READ(LUIN,*)DR4J_LUM(1:ND)
	      L_DR4J(I)=SUM(DR4J_LUM(1:ND))
	      DONL_DR4J_LUM=.TRUE.
	    END IF
	    IF(INDEX(STRING,'Energy deposited') .NE. 0 .AND. .NOT. DONL_DECAY_LUM)THEN
	      READ(LUIN,*)DECAY_LUM(1:ND)
	      L_DECAY(I)=SUM(DECAY_LUM(1:ND))
	      DONL_DECAY_LUM=.TRUE.
	    END IF
	  END DO
200	  CLOSE(UNIT=10)
!
! Get the SN's age from the VADAT file.
!
	  FILE_NAME=TRIM(DIR_NAME(I))//'VADAT'
	  OPEN(UNIT=LUIN,FILE=FILE_NAME,STATUS='OLD',ACTION='READ')
	  DO WHILE(SN_AGE(I) .EQ. 0)
	    READ(LUIN,'(A)')STRING
	    IF(INDEX(STRING,'[SN_AGE]') .NE. 0)THEN
	      READ(STRING,*)SN_AGE(I)
	      CLOSE(UNIT=10)
	    END IF
	  END DO
!
! Read J so that we can compute the raditive energy at each time step.
!
	  FILE_NAME=TRIM(DIR_NAME(I))//'JH_AT_CURRENT_TIME'
	  CALL READ_DIRECT_INFO_V3(J,REC_LENGTH,FILE_DATE,FILE_NAME,LUIN,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Error opening/reading INFO file: check format'
	    WRITE(6,*)'Also check error file or fort.2'
	  END IF
          OPEN(UNIT=LUIN,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',
	1         RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS)
	   READ(LUIN,REC=3)ST_REC,NCF,ND
	   ST_REC=ST_REC+2
	   IF(ALLOCATED(RJ))DEALLOCATE(R,RJ,HFLUX,NU)
	   ALLOCATE(RJ(ND,NCF))
	   ALLOCATE(HFLUX(ND,NCF))
	   ALLOCATE(NU(NCF))
	   ALLOCATE(R(ND))
!	   READ(LUIN,REC=ST_REC-2,IOSTAT=IOS)(R(1:ND))
	   READ(LUIN,REC=ST_REC-2)R(1:ND)
	   DO ML=1,NCF
	     READ(LUIN,REC=ST_REC+ML-1)RJ(1:ND,ML),HFLUX(1:ND-1,ML),T1,T2,NU(ML)
	   END DO
	   NCF=NCF-1
	   CLOSE(UNIT=10)
!
! Compare R in RVTJ with that in the JH file.
!
	  DO J=1,ND
	    IF( ABS(1.0D0-R_RVTJ(J)/R(J)) .GT. 1.0D-07)THEN
	      WRITE(6,*)'Error R grid does not match'
	      WRITE(6,*)R(J)
	      WRITE(6,*)R_RVTJ(J)
	      WRITE(6,*)'Directory name is ',DIR_NAME(I)
	      STOP
	    END IF
	  END DO
!
	  YV=0.0D0
          DO ML=1,NCF-1
            DO J=1,ND
              YV(J)=YV(J)+(NU(ML)-NU(ML+1))*(RJ(J,ML)+RJ(J,ML+1))
            END DO
          END DO
!
! We only print the last 15 characters of the directory name.
!
	  K=LEN_TRIM(DIR_NAME(I)); J=MAX(1,K-15)
	  WRITE(6,'(1X,A,T21,A,F9.4,3X,A,I4)')DIR_NAME(I)(J:K),'AGE(d)=',SN_AGE(I),'ND=',ND            !YV(1),YV(ND)
	  CALL LUM_FROM_ETA(YV,R,ND)
	  E_RAD(I)=2.0D0*PI*4.0D+37*PI*SUM(YV(1:ND))/C_MMS
!
	  WRITE(LUOUT,'(F10.4,7ES14.4,4X,A)')SN_AGE(I),L_OBS(I),L_CMF(I),L_MECH(I),L_DR4J(I),
	1               L_INTERN(I),L_DECAY(I),E_RAD(I),DIR_NAME(I)(MAX(1,K-10):K)
	  FLUSH(UNIT=LUOUT)
	END DO
!
	WRITE(6,'(A)')RED_PEN
	WRITE(6,*)'Output is to TIME_INTEGRATED_ENERGY_CHK'
	WRITE(6,'(A)')DEF_PEN
!
	SN_AGE(1:NMOD)=SN_AGE(1:NMOD)*24.0D0*3600.0D0
	T1=0.0D0; T2=0.0D0; T3=0.0D0
	WRITE(LUOUT,*)' '
	WRITE(LUOUT,*)' '
	WRITE(LUOUT,*)' The 6th and 7th (i.e., second last) columns should be identical'
	WRITE(LUOUT,*)' The % Change is the increase in Int(L-Q+I)t+t.E over a single time step.'
	WRITE(LUOUT,*)' The % Error  is the increase in Int(L-Q+I)t+t.E over from the first time step.'
	WRITE(LUOUT,*)' The % E2     is 200*dE/ABS(all terms) -  gives depature from energy conservation irrespective of dominant E'
	WRITE(LUOUT,*)' '
	WRITE(LUOUT,'(4X,A,4(4X,A),X,A,5X,A,3(4X,A))')'AGE(d)','Int(tL.dt)','Int(tQ.dt)','Int(tI.dt)','  t.E(rad)',
	1              'Int(L-Q+I)t+t.E','t(1).E(1)','% Change', ' % Error',' % E'
	E_CONS(1)=SN_AGE(1)*E_RAD(1)
	DO I=2,NMOD
	   T1=T1+(SN_AGE(I)-SN_AGE(I-1))*(SN_AGE(I)*L_CMF(I)+SN_AGE(I-1)*L_CMF(I-1))*0.5D0*3.826D+33
	   T2=T2+(SN_AGE(I)-SN_AGE(I-1))*(SN_AGE(I)*L_DECAY(I)+SN_AGE(I-1)*L_DECAY(I-1))*0.5D0*3.826D+33
	   T3=T3+(SN_AGE(I)-SN_AGE(I-1))*(SN_AGE(I)*L_INTERN(I)+SN_AGE(I-1)*L_INTERN(I-1))*0.5D0*3.826D+33
	   E_CONS(I)=T1-T2+T3+SN_AGE(I)*E_RAD(I)
	   T5=100.0D0*(E_CONS(I)-SN_AGE(1)*E_RAD(1))
	   T4=T5/SN_AGE(1)/E_RAD(1)
	   T5=2.0D0*T5/(ABS(T1)+ABS(T2)+ABS(T3)+SN_AGE(I)*E_RAD(I)+SN_AGE(1)*E_RAD(1))
	   K=LEN_TRIM(DIR_NAME(I)); J=MAX(1,K-15)
	   WRITE(LUOUT,'(F10.4,4ES14.4,ES16.4,ES14.4,3ES12.2,3X,A)')SN_AGE(I)/24.0D0/3600.0D0,
	1              T1,T2,T3,SN_AGE(I)*E_RAD(I),
	1              E_CONS(I),SN_AGE(1)*E_RAD(1),
	1              100.0D0*(E_CONS(I)-E_CONS(I-1))/SN_AGE(1)/E_RAD(1),
	1              T4,T5,DIR_NAME(I)(MAX(1,K-10):K)
	END DO
!
! Same quantities as above, but a different grouping.
!
	T1=0.0D0; T2=0.0D0; T3=0.0D0
	WRITE(LUOUT,*)' '
	WRITE(LUOUT,*)' '
	WRITE(LUOUT,*)' The last two columns should be compared.'
	WRITE(LUOUT,*)' '
	WRITE(LUOUT,'(78X,A)')'t(1).E(1)+'
	WRITE(LUOUT,'(4X,A,4(3X,A),A,4X,A)')'AGE(d)','Int(tL.dt)','Int(tQ.dt)','Int(tI.dt)','  t.E(rad)',
	1              '  Int(tL)+t.E','Int(Q-I)t'
	E_CONS(1)=SN_AGE(1)*E_RAD(1)
	DO I=2,NMOD
	   T1=T1+(SN_AGE(I)-SN_AGE(I-1))*(SN_AGE(I)*L_CMF(I)+SN_AGE(I-1)*L_CMF(I-1))*0.5D0*3.826D+33
	   T2=T2+(SN_AGE(I)-SN_AGE(I-1))*(SN_AGE(I)*L_DECAY(I)+SN_AGE(I-1)*L_DECAY(I-1))*0.5D0*3.826D+33
	   T3=T3+(SN_AGE(I)-SN_AGE(I-1))*(SN_AGE(I)*L_INTERN(I)+SN_AGE(I-1)*L_INTERN(I-1))*0.5D0*3.826D+33
	   E_CONS(I)=T1+SN_AGE(I)*E_RAD(I)
	   WRITE(LUOUT,'(F10.4,6ES13.4)')SN_AGE(I)/24.0D0/3600.0D0,
	1              T1,T2,T3,SN_AGE(I)*E_RAD(I),
	1              E_CONS(I),SN_AGE(1)*E_RAD(1)+T2-T3
	END DO
!
! We directly integrate equation (50) in Hillier and Dessart (2012) over time.
! This gives an alterative energy constraint that is useful for checking the modeling, but is less
! useful observationlly.
!
	T1=0.0D0; T2=0.0D0; T3=0.0D0
	WRITE(LUOUT,*)' '
	WRITE(LUOUT,*)' '
	WRITE(LUOUT,*)' The last columns should be compared.'
	WRITE(LUOUT,*)' Alternatively compare Int(L.dt) with (Q-I+S+E1-EI).'
	WRITE(LUOUT,*)' '
	WRITE(LUOUT,'(79X,A)')' S='
	WRITE(LUOUT,'(4X,A,5(4X,A),2X,A,2X,A,X,A,5X,A)')
	1              'AGE(d)','Int(L.dt)','Int(Q.dt)','Int(I.dt)','E(rad)(1)',
	1              'E(rad)(I)','Int(E/t.dt)','Q-I-S+E1-EI','Int(L.dt)+EI','Q-I+S+EI'
	E_CONS(1)=0.0D0
	DO I=2,NMOD
	   T1=T1+(SN_AGE(I)-SN_AGE(I-1))*(L_CMF(I)+L_CMF(I-1))*0.5D0*3.826D+33
	   T2=T2+(SN_AGE(I)-SN_AGE(I-1))*(L_DECAY(I)+L_DECAY(I-1))*0.5D0*3.826D+33
	   T3=T3+(SN_AGE(I)-SN_AGE(I-1))*(L_INTERN(I)+L_INTERN(I-1))*0.5D0*3.826D+33
	   E_CONS(I)=E_CONS(I-1)+(SN_AGE(I)-SN_AGE(I-1))*(E_RAD(I)/SN_AGE(I)+E_RAD(I-1)/SN_AGE(I-1))*0.5D0
	   WRITE(LUOUT,'(F10.4,9ES13.4)')SN_AGE(I)/24.0D0/3600.0D0,
	1              T1,T2,T3,E_RAD(1),E_RAD(I),
	1              E_CONS(I),T2-T3+(E_RAD(1)-E_RAD(I)-E_CONS(I)),
	1              T1+E_RAD(I),T2-T3+(E_RAD(1)-E_CONS(I))
	END DO
!
	STOP
	END
