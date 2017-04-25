!
! Routine to pack levels: Packing can be done for
!          States with LS nomenclature.
!          States with JL nomenclature.
!
! Needs oscillator file
!
	PROGRAM PACK_OSC
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER NLEV
	INTEGER IOS
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: LUIN=7
	INTEGER, PARAMETER :: IZERO=0
	INTEGER I
!
	CHARACTER(LEN=10) ION_ID
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=132) OSC_FILE
!
! Cleaned 23-Dec-2004
! Created 02-Jun-2003
!
	WRITE(T_OUT,'()')
	WRITE(T_OUT,'(70A)')('*',I=1,70)
	WRITE(T_OUT,'()')
	WRITE(T_OUT,'(A)')' The name of the oscillator file is prompted.'
	WRITE(T_OUT,'()')
	WRITE(T_OUT,'(A)')' The following diagnostic file are output:'
	WRITE(T_OUT,'(A)')'                  XzV_PACK,'
	WRITE(T_OUT,'(A)')'          SLP_CHK_FOR_XzV,'
	WRITE(T_OUT,'(A)')'         ERROR_CHK_FOR_XzV,'
	WRITE(T_OUT,'()')
	WRITE(T_OUT,'(70A)')('*',I=1,70)
	WRITE(T_OUT,'()')
!
	ION_ID=' '
	CALL GEN_IN(ION_ID,'Ionization identification (e.g., CIV)')
!
! Get number of levels in Oscilator file.
!
	IOS=1
	DO WHILE(IOS .NE. 0)
	  OSC_FILE=' '
	  CALL GEN_IN(OSC_FILE,'File with oscilator strengths')
	  CALL GEN_ASCI_OPEN(LUIN,OSC_FILE,'OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error occurred reading Oscillator file: try again'
	  END IF
	END DO
	STRING=' '
	DO WHILE( INDEX(STRING,'!Number of energy levels') .EQ. 0)
	  READ(LUIN,'(A)')STRING
	END DO
	CLOSE(LUIN)
	READ(STRING,*)NLEV
	CLOSE(LUIN)
	WRITE(T_OUT,*)'Number of atomic levels is',NLEV
!
	CALL PACK_OSC_SUB(OSC_FILE,ION_ID,NLEV)
!
	STOP
	END
!
!
!
	SUBROUTINE PACK_OSC_SUB(OSC_FILE,ION_ID,NLEV)
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER NLEV
	CHARACTER(LEN=*) ION_ID
	CHARACTER(LEN=*) OSC_FILE
!
! Atomic data variables for main species.
!
	REAL*8 FEDGE(NLEV)
	REAL*8 ENERGY(NLEV)
	REAL*8 G(NLEV)
	REAL*8 LAM_EDGE(NLEV)
	REAL*8 E_STRT(NLEV)
!
	REAL*8 EDGE_SUM(NLEV)
	REAL*8 G_SUM(NLEV)
	REAL*8 FOSC(NLEV,NLEV)
!
	REAL*8 GF_CUT
	INTEGER GF_LEV_CUT
	INTEGER MIN_NUM_TRANS
!
	INTEGER CROSS_TYPE(NLEV)
	INTEGER F_TO_S(NLEV)
	CHARACTER*30 NAME(NLEV)
	CHARACTER*80 FIN_STATE
	CHARACTER*80 FILENAME
!
	REAL*8 AT_NO
	REAL*8 ZION
	REAL*8 ZION_RD
	REAL*8 GION
	REAL*8 PHOT_GION
	REAL*8 EXC_EN
	REAL*8 ION_EN
	CHARACTER*20 EN_DATE
!
	LOGICAL AF_CHK
	LOGICAL PACK
	LOGICAL USE_G_WEIGHTING
!
! Variables to read in dielectronic lines.
!
	LOGICAL DO_DIE
	LOGICAL DO_DIE_REG
	LOGICAL DO_DIE_WI
!
! Used for dynamic smoothing of th ephotoioization cross-sections.
!
	REAL*8 VSM_DIE_KMS
	REAL*8 SIG_GAU_KMS
	REAL*8 FRAC_SIG_GAU
	REAL*8 CUT_ACCURACY
	LOGICAL ABOVE_EDGE
!
	CHARACTER*30 LS_NAME(NLEV)	!Term (LS) designation
	CHARACTER*1 LEV_ANG(NLEV)	!Total angular momentum
	CHARACTER*1 LEV_SPIN(NLEV)	!Multiplicity
	CHARACTER*1 LEV_PARITY(NLEV)	!Parity
!
! Functions to return parity and multiplicity..
!
	CHARACTER*1 PARITY
	CHARACTER*1 SPIN
!
	LOGICAL PACK_LJ_STATES
	REAL*8 EMIN_PACK
	INTEGER IMIN_PACK
!
! For atom packed according to terms (generally LS).
! We allocate NLEV for NAME_PACK, since NLEV is a maximum, 
! an we ned the storage before we know how many packed levels
! we have.
!
	CHARACTER(LEN=30) NAME_PACK(NLEV)
	REAL*8, ALLOCATABLE :: EDGE_PACK(:)
	REAL*8, ALLOCATABLE :: FOSC_PACK(:,:)
	REAL*8, ALLOCATABLE :: G_PACK(:)
	REAL*8, ALLOCATABLE :: RECOM_PACK(:,:)
	REAL*8, ALLOCATABLE :: ARAD(:)
	REAL*8, ALLOCATABLE :: GAM2(:)
	REAL*8, ALLOCATABLE :: GAM4(:)
	LOGICAL, ALLOCATABLE :: KNOWN_ENERGY_LEVEL(:)
	LOGICAL, ALLOCATABLE :: SECND(:,:)
	INTEGER, ALLOCATABLE :: TRANS(:,:)
	INTEGER, ALLOCATABLE :: CROSS_TYPE_PACK(:)
!
	REAL*8, ALLOCATABLE :: FLOW_SUM(:)
	REAL*8, ALLOCATABLE :: FHIGH_SUM(:)
	REAL*8, ALLOCATABLE :: ALOW_SUM(:)
	REAL*8, ALLOCATABLE :: AHIGH_SUM(:)
	LOGICAL, ALLOCATABLE :: LEVEL_DONE(:)
!
	REAL*8, ALLOCATABLE :: FOSC_TMP(:,:)
	REAL*8, ALLOCATABLE :: VEC_DP_WRK(:)
	INTEGER, ALLOCATABLE :: VEC_INDX(:)
!
	INTEGER, PARAMETER :: NANG=12
	CHARACTER*1 ANG_STR(NANG)
	DATA ANG_STR/'S','P','D','F','G','H','I','K','L','M','N','O'/         !,'Z','W'/
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN               
	DOUBLE PRECISION CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	CHARACTER(LEN=1) FORMFEED
	LOGICAL FILE_OPEN
!
	INTEGER IZERO
	PARAMETER (IZERO=0)
!
	INTEGER I,J,K,L,IS,JS,IOS
	INTEGER ID,PHOT_ID
	INTEGER CNT
	INTEGER N_LS_TERMS
	INTEGER MAX_NAME_LNGTH
	REAL*8 T1,T2
!
	CHARACTER(LEN=30) TMP_STR
	CHARACTER(LEN=24) TIME
	CHARACTER(LEN=80) STRING
!
	LOGICAL L_TRUE,L_FALSE
	DATA L_TRUE/.TRUE./
	DATA L_FALSE/.FALSE./
!
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6		!For terminal IO              
!
	INTEGER ERROR_LU
	INTEGER LUER
	INTEGER, PARAMETER :: LUIN=30			!File IO
	INTEGER, PARAMETER :: LUOUT=40
!
	EXTERNAL SPIN
	EXTERNAL PARITY
	EXTERNAL ERROR_LU
!
! 
! Set constants.
!
	CHIBF=2.815D-06
	CHIFF=3.69D-29
	HDKT=4.7994145D0
	TWOHCSQ=0.0147452575D0
	OPLIN=2.6540081D+08
	EMLIN=5.27296D-03
	FORMFEED=''
!
	LUER=2      					!ERROR_LU()
	FILENAME='ERROR_CHK_FOR_'//TRIM(ION_ID)
	OPEN(UNIT=LUER,FILE=FILENAME,STATUS='UNKNOWN')
!
! Read in all the atomic data. Some of these are superfluous.
!
	GF_CUT=0.0D0			!These ensure we get all transitions.
	GF_LEV_CUT=NLEV+1
	MIN_NUM_TRANS=NLEV*NLEV
	EN_DATE=' '
        CALL RD_ENERGY(NAME,G,ENERGY,FEDGE,NLEV,NLEV,			!To get ENERGY
	1       ION_EN,ZION,EN_DATE,OSC_FILE,LUIN,LUOUT,IOS)
	CLOSE(LUIN)
	CALL GENOSC_V6(FOSC,FEDGE,G,NAME,ION_EN,ZION,EN_DATE,NLEV,I,
	1        'SET_ZERO',GF_CUT,GF_LEV_CUT,MIN_NUM_TRANS,
	1        LUIN,LUOUT,OSC_FILE)
!
!
! Check that the states with [] have correct statistical weight.
!
	WRITE(LUER,*)' '
	DO I=1,NLEV
	  MAX_NAME_LNGTH=MAX(MAX_NAME_LNGTH,LEN_TRIM(NAME(I)))
	  J=INDEX(NAME(I),'[')
	  IF(J .EQ. 0)THEN
	  ELSE
	    K=INDEX(NAME(I),']')
	    TMP_STR=NAME(I)(J+1:K-1)
	    J=INDEX(TMP_STR,'/2')
	    IF(J .EQ. 0)THEN
	      READ(TMP_STR,*)K
	      K=2*K+1
	    ELSE
	      READ(TMP_STR(1:J-1),*)K
	      K=K+1
	    END IF
	    IF(K .NE. NINT(G(I)))THEN
	       WRITE(T_OUT,*)'Error --- invalid statistical weight for level ',NAME(I)
	       WRITE(T_OUT,*)'Statistical weight is',G(I)
	       WRITE(T_OUT,*)'Computed weight is',K
	       WRITE(LUER,*)'Error --- invalid statistical weight for level ',NAME(I)
	       WRITE(LUER,*)'Statistical weight is',G(I)
	       WRITE(LUER,*)'Computed weight is',K
	    END IF
	  END IF
	END DO
!
! Strip J values from names, and determine parity and spin.
!
	MAX_NAME_LNGTH=0
	DO I=1,NLEV
	  MAX_NAME_LNGTH=MAX(MAX_NAME_LNGTH,LEN_TRIM(NAME(I)))
	  J=INDEX(NAME(I),'[')
	  IF(J .EQ. 0)THEN
	    LS_NAME(I)=NAME(I)
	  ELSE
	    LS_NAME(I)=NAME(I)(1:J-1)
	  END IF
	  LEV_PARITY(I)=PARITY(LS_NAME(I))
	  LEV_SPIN(I)=SPIN(LS_NAME(I))
	END DO
!
! Determine number of LS terms in FULL model atom.
!
	N_LS_TERMS=1
	DO I=2,NLEV
	  J=1
	  DO WHILE(LS_NAME(I) .NE. LS_NAME(J) .AND. J .LE. I-1)
	    J=J+1
	  END DO
	  IF(J .EQ. I)N_LS_TERMS=N_LS_TERMS+1
	END DO
!
! Check that states have correct statistical weight.
!
	WRITE(LUER,*)' '
	DO I=1,NLEV
	 LEV_ANG(I)=' '
	  K=LEN_TRIM(NAME(I))
	  J=INDEX(NAME(I),'[')
	  IF(J .NE. 0)K=J-1
	  IF(NAME(I)(K-2:K) .EQ. 'SNG')THEN
	  ELSE IF(NAME(I)(K-2:K) .EQ. 'TRP')THEN
	  ELSE IF(NAME(I)(K-2:K) .EQ. '___')THEN
	  ELSE
	    K=K-1
	    IF(NAME(I)(K:K) .GE. 'a' .AND. NAME(I)(K:K) .LE. 'z')K=K-1
	    IF(NAME(I)(K:K) .EQ. '1' .OR. NAME(I)(K:K) .LE. '2')THEN
	      K=K-1
	    END IF
	    IF(NAME(I)(K:K) .EQ. 'W' .OR. NAME(I)(K:K) .EQ. 'Z')THEN
	      LEV_ANG(I)=NAME(I)(K:K)
	    ELSE
	      T1=0.0D0
	      DO J=1,NANG
	        IF(NAME(I)(K:K) .EQ. ANG_STR(J))THEN
	          LEV_ANG(I)=ANG_STR(J)
	          IF(INDEX(NAME(I),'[') .EQ. 0)THEN
	            READ(LEV_SPIN(I),*)T1
	            T1=T1*(2*J-1)
	            IF(NINT(T1) .NE. NINT(G(I)))THEN
	              WRITE(LUER,*)' '
	              WRITE(LUER,*)'Invalid statistical weight for state',NAME(I)
	              WRITE(T_OUT,*)'Invalid statistical weight for state',NAME(I)
	            END IF
	          END IF
	          EXIT
	        END IF
	      END DO
	      IF(LEV_ANG(I) .EQ. ' ')THEN
	         WRITE(LUER,*)'Unable to get L for ',TRIM(NAME(I))
	         WRITE(T_OUT,*)'Unable to get L for ',TRIM(NAME(I))
	      END IF
	    END IF
	  END IF
	END DO
!
! Create a summary file with f and A values summed over transitions to individual leveles
!
	AF_CHK=.FALSE.
	CALL GEN_IN(AF_CHK,'Create check file with A & F summed over levels?')
	IF(AF_CHK)THEN
	  FILENAME='AF_CHK_FOR_'//TRIM(ION_ID)
	  OPEN(UNIT=30,STATUS='UNKNOWN',FILE=FILENAME)
	  WRITE(30,'()')
	  WRITE(30,'()')
	  WRITE(30,'(A)')' Summary file with information on oscillator strength'//
	1                     ' and Einstein A values.'
	  WRITE(30,'(A)')'   FL_SUM is the sum of f from lower state to the given state.'
	  WRITE(30,'(A)')'   AL_SUM is the sum of the decay rates from the given state.'
	  WRITE(30,'(A)')'   FH_SUM is the sum of f from higher states to the given state.'
	  WRITE(30,'(A)')'   AH_SUM is the sum of A from higher states to the given state.'
	  WRITE(30,'(A)')' These sums may be weighted by the statistical weights.'
	  USE_G_WEIGHTING=.FALSE.
	  CALL GEN_IN(USE_G_WEIGHTING,'Weight f,A sums by g?')
	  WRITE(30,'()')
	  STRING='Name'
	  
	  IF(USE_G_WEIGHTING)THEN
	    WRITE(30,'(A,3X,3(1X,A),8X,A,3X,4(3X,A,4X))')STRING(1:MAX_NAME_LNGTH),
	1           'S','L','P','G','GFL_SUM','GFH_SUM','GAL_SUM','GAH_SUM'
	  ELSE
	     WRITE(30,'(A,3X,3(1X,A),8X,A,3X,4(4X,A,4X))')STRING(1:MAX_NAME_LNGTH),
	1           'S','L','P','G','FL_SUM','FH_SUM','AL_SUM','AH_SUM'
	  END IF
	  WRITE(30,'()')
!
	  ALLOCATE (FLOW_SUM(NLEV))
	  ALLOCATE (ALOW_SUM(NLEV))
	  ALLOCATE (FHIGH_SUM(NLEV))
	  ALLOCATE (AHIGH_SUM(NLEV))
	  FLOW_SUM=0; FHIGH_SUM=0
	  ALOW_SUM=0; AHIGH_SUM=0
!
	  DO I=1,NLEV
	    IF(USE_G_WEIGHTING)THEN
	      FLOW_SUM(I)=SUM(G(1:I-1)*FOSC(1:I-1,I))
	      ALOW_SUM(I)=SUM(G(I)*FOSC(I,1:I-1))
	      FHIGH_SUM(I)=SUM(G(I)*FOSC(I,I+1:NLEV))
	      AHIGH_SUM(I)=SUM(G(I+1:NLEV)*FOSC(I+1:NLEV,I))
	    ELSE
	      FLOW_SUM(I)=SUM(FOSC(1:I-1,I))
	      ALOW_SUM(I)=SUM(FOSC(I,1:I-1))
	      FHIGH_SUM(I)=SUM(FOSC(I,I+1:NLEV))
	      AHIGH_SUM(I)=SUM(FOSC(I+1:NLEV,I))
	    END IF
	    WRITE(30,'(A,3X,3(A,2X),3X,F5.0,4ES14.4)')NAME(I)(1:MAX_NAME_LNGTH),
	1                  LEV_SPIN(I),LEV_ANG(I),LEV_PARITY(I),
	1                  G(I),FLOW_SUM(I),FHIGH_SUM(I),ALOW_SUM(I),AHIGH_SUM(I)
	  END DO
	  CLOSE(UNIT=30)
	END IF
!
	PACK=.TRUE.
	PACK_LJ_STATES=.TRUE.
	EMIN_PACK=-1.0D0			!i.e. pack all states
	CALL GEN_IN(PACK,'Pack states?')
	CALL GEN_IN(EMIN_PACK,'Pack when states above this energy (in cm^{-1})')
	CALL GEN_IN(PACK_LJ_STATES,'Pack states in intermediate coupling?')
!
	ALLOCATE (LEVEL_DONE(NLEV))
	LEVEL_DONE=.FALSE.
!
! Write summary of levels with S, L & parity. Levels are groupe according to
! their multiplet name.
!
	IF(PACK)THEN
	  FILENAME='SLP_CHK_FOR_'//TRIM(ION_ID)
	  OPEN(UNIT=30,STATUS='UNKNOWN',FILE=FILENAME)
	  WRITE(30,'()')
	  STRING='Level'
	  WRITE(30,'(A,3X,3(2X,A))')STRING(1:MAX_NAME_LNGTH),'S','L','P'
	  DO J=1,NLEV
	    IF(.NOT. LEVEL_DONE(J))THEN
	      DO I=J,NLEV
	        IF(LS_NAME(J) .EQ. LS_NAME(I))THEN
	          WRITE(30,'(A,3X,3(2X,A))')NAME(I)(1:MAX_NAME_LNGTH),
	1                                LEV_SPIN(I),LEV_ANG(I),LEV_PARITY(I)
	          LEVEL_DONE(I)=.TRUE.
	        END IF
	      END DO
	    END IF
	  END DO
	END IF
!
! Pack levels
!
	IF(PACK)THEN
	  F_TO_S(:)=0
	  CNT=1
!
	  DO I=1,NLEV
	    WRITE(41,*)I,ENERGY(I),EMIN_PACK
	    IF(ENERGY(I) .GT. EMIN_PACK)THEN
	      IMIN_PACK=I
	      EXIT
	    END IF
	  END DO
	  WRITE(6,'(A,I5,A)')'Packing levels',IMIN_PACK,' and above.'
!
	  CNT=0
	  DO I=1,IMIN_PACK-1
	    NAME_PACK(I)=NAME(I)
	    F_TO_S(I)=I
	    CNT=I
	  END DO
!
	  IF(PACK_LJ_STATES)THEN
	    DO I=IMIN_PACK,NLEV
	      K=INDEX(LS_NAME(I),'{')
	      IF(K .NE. 0)THEN
	        J=INDEX(LS_NAME(I)(K:),'e')
	        IF(J .EQ. 0)THEN
	          LS_NAME(I)=LS_NAME(I)(1:K)//'}o'
	        ELSE
	          LS_NAME(I)=LS_NAME(I)(1:K)//'}e'
	        END IF
	      END IF
	    END DO
	  END IF
!
	  WRITE(T_OUT,'(A)')' Will now pack the levels'
	  DO I=IMIN_PACK,NLEV
	    J=INDEX(NAME(I),'[')
	    K=INDEX(NAME(I),'{')
	    IF(J .EQ. 0 .AND. K .EQ. 0)THEN
	      CNT=CNT+1
	      F_TO_S(I)=CNT
	      NAME_PACK(CNT)=LS_NAME(I)
	    ELSE
	      J=1
	      DO WHILE(F_TO_S(I) .EQ. 0 .AND. J .LE. I-1)
	        IF( LS_NAME(I) .EQ. LS_NAME(J) )THEN
	          F_TO_S(I)=F_TO_S(J)
	          EXIT
	        END IF
	        J=J+1
	      END DO
	      IF(F_TO_S(I) .EQ. 0)THEN
	        CNT=CNT+1
	        F_TO_S(I)=CNT
	        NAME_PACK(CNT)=LS_NAME(I)
	      END IF
	    END IF
	  END DO
!
! Can now do packing.
!
	  ALLOCATE (EDGE_PACK(CNT))
	  ALLOCATE (G_PACK(CNT))
	  ALLOCATE (FOSC_PACK(CNT,CNT))
	  EDGE_PACK(:)=0.0D0; G_PACK(:)=0.0D0; FOSC_PACK(:,:)=0.0D0
	  DO I=1,NLEV
	    J=F_TO_S(I)
	    G_PACK(J)=G_PACK(J)+G(I)
	    EDGE_PACK(J)=EDGE_PACK(J)+G(I)*FEDGE(I)
	  END DO
	  DO J=1,CNT
	    EDGE_PACK(J)=EDGE_PACK(J)/G_PACK(J)
	  END DO	
	  DO J=2,NLEV
	    JS=F_TO_S(J)
	    DO I=1,J
	      IS=F_TO_S(I)
	      FOSC_PACK(IS,JS)=FOSC_PACK(IS,JS)+G(I)*FOSC(I,J)
	      FOSC_PACK(JS,IS)=FOSC_PACK(JS,IS)+G(J)*FOSC(J,I)
	    END DO
	  END DO
	  DO JS=2,CNT
	    DO IS=1,JS-1
	    FOSC_PACK(IS,JS)=FOSC_PACK(IS,JS)/G_PACK(IS)
	    FOSC_PACK(JS,IS)=FOSC_PACK(JS,IS)/G_PACK(JS)
	   END DO
	  END DO
!
! Since we have packed the levels, some levels may be out of order.
! Therfore need to do a sort.
!
	  ALLOCATE (VEC_INDX(CNT))
	  ALLOCATE (VEC_DP_WRK(CNT))
	  ALLOCATE (FOSC_TMP(CNT,CNT))
!
	  WRITE(6,*)'Sorting levels into energy order'
	  CALL INDEXX(CNT,EDGE_PACK,VEC_INDX,L_FALSE)
	  CALL SORTDP(CNT,EDGE_PACK,VEC_INDX,VEC_DP_WRK)
	  CALL SORTDP(CNT,G_PACK,VEC_INDX,VEC_DP_WRK)
	  CALL SORTCHAR(CNT,NAME_PACK,VEC_INDX,NAME)
	  FOSC_TMP=FOSC_PACK
	  DO I=1,CNT
	    FOSC_PACK(I,:)=FOSC_TMP(VEC_INDX(I),:)
	  END DO
	  FOSC_TMP=FOSC_PACK
	  DO I=1,CNT
	    FOSC_PACK(:,I)=FOSC_TMP(:,VEC_INDX(I))
	  END DO
	  WRITE(6,*)'Finished sorting levels'
!
	DO I=1,CNT
	  WRITE(77,'(A)')TRIM(NAME_PACK(I))
	END DO
!
	  ALLOCATE (KNOWN_ENERGY_LEVEL(CNT))
	  ALLOCATE (ARAD(CNT))
	  ALLOCATE (GAM2(CNT))
	  ALLOCATE (GAM4(CNT))
	  ALLOCATE (TRANS(CNT,CNT))
	  ALLOCATE (SECND(CNT,CNT))
	  ARAD=0.0D0; GAM2=0.0D0; GAM4=0.0D0
	  KNOWN_ENERGY_LEVEL=.TRUE.
	  FILENAME=TRIM(ION_ID)//'_PACK'
	  CALL WRITE_OSC_V2(FOSC_PACK,TRANS,SECND,EDGE_PACK,G_PACK,
	1            NAME_PACK,ARAD,GAM2,GAM4,KNOWN_ENERGY_LEVEL,
	1            ION_EN,ZION,LUIN,CNT,L_FALSE,L_TRUE,EN_DATE,
	1            ION_ID,FILENAME)
!
	END IF			!End pack
!
	STOP
	END
!
!
	FUNCTION PARITY(NAME)
	IMPLICIT NONE
!
	CHARACTER*(*) NAME
	CHARACTER*1 PARITY
!
	EXTERNAL ICHRLEN
	INTEGER ICHRLEN,J,T_OUT
!
	T_OUT=6
	J=INDEX(NAME,'[')-1
	IF(J .LE. 0)J=ICHRLEN(NAME)
	PARITY=NAME(J:J)
	IF(PARITY .NE. 'e' .AND. PARITY .NE. 'o')THEN
	   IF(PARITY .NE. 'Z' .AND. PARITY .NE. '_')THEN
	     WRITE(T_OUT,*)'Warning: Indeterminate parity for level',NAME
	   END IF
	   PARITY=' '
	END IF
	RETURN
	END
!
!
	FUNCTION SPIN(NAME)
	IMPLICIT NONE
!
! Altered 31-Jan-2003: T_OUT set to 6 (instead of 5)
! Altered 20-Oct-2000: Spin determination altered to allow for names
!                        like .. 3Pbe or 3P2e.
!
	CHARACTER*(*) NAME
	CHARACTER*1 SPIN
C
	EXTERNAL ICHRLEN
	INTEGER ICHRLEN,LS,J,IOS,T_OUT
	T_OUT=6
C
	J=INDEX(NAME,'[')-1
	IF(J .LE. 0)J=ICHRLEN(NAME)
	IF(NAME(J:J) .EQ. 'e' .OR. NAME(J:J) .eq. 'o')THEN
	  LS=J-2
	ELSE
	  LS=J-1
	END IF
!
! Allow for names of the form ... 3Pbe or 3P2e where b and 2 are used to
! break name degeneracies.
!
	IF(NAME(LS:LS) .GT. 'A' .AND. NAME(LS:LS) .LE. 'Z')LS=LS-1
!
	SPIN=NAME(LS:LS)
	IOS=0
	READ(SPIN,'(I1)',IOSTAT=IOS)J
	IF(IOS .NE. 0)THEN
	   WRITE(T_OUT,*)'Indeterminate spin for level',NAME
	   SPIN=' '
	END IF
	RETURN
	END	
