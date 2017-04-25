!
! General routine to return collsion rates (in the form of collsion strengths
! OMEGA) for a single depth in the atmosphere. Two methods are used:
!	(1) Values are read in from a file. This file will not be read
!                in for every depth UNLESS it has been corrupted.
!                These values are used unless unavailable.
!                Interpolation is done in the Log-Log plane.
!	 (2) Use the approximate formulae which express the collsion
!                strengths interms of the oscilator values.
!
	SUBROUTINE OMEGA_GEN_V3(OMEGA,dln_OMEGA_dlnT,EDGE,EIN_A,STAT_WT,
	1                       LEVELNAME,ZION,NLEV,TEMP,
	1                       ID,FILE_NAME)
	IMPLICIT NONE
!
! Altered 04-Oct-2016 : Efectively we now only compute the threshold photoionzation cross sections
!                         on the first entry (change was done in test routine earlier). 
! Altered 12-Oct-2012 : MAX_TRANS increased to 50,000
! Altered 23-Nov-2007 : MAX_TVALS adjusted upwards from 21 to 40 (20-Nov-2007).
! Altered 23-Feb-1999 : Access of collison table increased to improve speed
!                         for Fe2. Old method very inefficient. Now loop
!                         over TABLE entries rather than NL(I) and NUP(J).
! Altered 08-Feb-1999 : MAX_TRANS increased from 500 to 2000 (for Fe2).
! Altered 05-Dec-1996 : PHOT_FUN replaced by PHOT_SUB. All photoioization
!                        cross sections computed at once. Called V2.
!                        STAT_WT include in call to GEN_OMEGA_RD(_V2)
! Altered 20-Jun-1996 : call to DP_ZERO removed.
!                       LST_RD_FILE included again. (16-Jun-1996)
! Altered 13-Jun-1996 : OMEGA_SET, OMEGA_SCALE now saved. OMEGA_SET is now
!                         taken from collisional file, and not put to 0.1
!                         always.
!
! Altered 25-May-1996 - THREE inserted.
! Created 13-Sep-1995
!
	INTEGER NLEV
	REAL*8 OMEGA(NLEV,NLEV),dln_OMEGA_dlnT(NLEV,NLEV)
!
	REAL*8 EIN_A(NLEV,NLEV)		!Einstein A coefficient
	REAL*8 EDGE(NLEV)		!Ionization frequency (10^15 Hz)
	REAL*8 STAT_WT(NLEV)		!Statistical weight.
	REAL*8 ZION
	CHARACTER*(*) LEVELNAME(NLEV)
!
	REAL*8 TEMP	       		!Temperature (10^4 K)
	CHARACTER*(*) FILE_NAME
	INTEGER ID                    !Indicates species for photoiozation data
!
	INTEGER ERROR_LU
	LOGICAL SAME_N
	EXTERNAL ERROR_LU
!
	REAL*8, PARAMETER :: THREE=3.0D0
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER MAX_TRANS,MAX_TVALS,MAX_TAB_SIZE
	PARAMETER (MAX_TRANS=50000)
	PARAMETER (MAX_TVALS=40)
	PARAMETER (MAX_TAB_SIZE=MAX_TVALS*MAX_TRANS)
	INTEGER, SAVE :: ID_LOW(MAX_TRANS)
	INTEGER, SAVE :: ID_UP(MAX_TRANS)
	INTEGER, SAVE :: ID_INDX(MAX_TRANS)
	REAL*8, SAVE :: T_TABLE(MAX_TVALS)
	REAL*8, SAVE :: OMEGA_TABLE(MAX_TAB_SIZE)
	REAL*8, SAVE :: OMEGA_SET,OMEGA_SCALE
	INTEGER, SAVE :: NUM_TRANS,NUM_TVALS
 
	CHARACTER*80 LAST_RD_FILE
	DATA LAST_RD_FILE/' '/
	SAVE LAST_RD_FILE
!
! Local variables.
!
	INTEGER I		!Used for lower level.
	INTEGER J		!Used for upper level.
	INTEGER PHOT_ID
	INTEGER L,K,LUER,IFAIL
!
	REAL*8 EX_E1X,EX_E1X_FUN
	REAL*8 X,FL,G1,G2,GBAR,dln_GBAR_dlnT,PHOT_CONST
	REAL*8 ALPHA
	REAL*8 T1
!
	REAL*8, PARAMETER :: ZERO=0.0D0
	REAL*8, SAVE, ALLOCATABLE :: PHOT_CROSS(:)
	LOGICAL RET_EDGE_CROSS
!
! Read in the atomic data. The data is not read in if the filename matches
! the filename of the previusly read data.
!
	CALL TUNE(1,'OMEGA_RD')
	IF(LAST_RD_FILE .NE. FILE_NAME)THEN
	  CALL GEN_OMEGA_RD_V2(OMEGA_TABLE,T_TABLE,
	1                      ID_LOW,ID_UP,ID_INDX,
	1                      OMEGA_SCALE,OMEGA_SET,
	1                      LEVELNAME,STAT_WT,NLEV,FILE_NAME,
	1                      NUM_TRANS,NUM_TVALS,
	1                      MAX_TRANS,MAX_TVALS,MAX_TAB_SIZE)
	END IF
	CALL TUNE(2,'OMEGA_RD')
	
!
! 
!
	LUER=ERROR_LU()
	OMEGA(:,:)=0.0D0			!NLEV*NLEV
	dln_OMEGA_dlnT(:,:)=0.0D0		!NLEV*NLEV
!
	IF(NUM_TRANS .EQ. 0)GOTO 1000
	IF(TEMP .LE. T_TABLE(1))THEN
	  DO K=1,NUM_TRANS
	    I=ID_LOW(K)
	    J=ID_UP(K)
	    OMEGA(I,J)=OMEGA_TABLE(ID_INDX(K))
	  END DO
	ELSE IF(TEMP .GE. T_TABLE(NUM_TVALS))THEN
	  DO K=1,NUM_TRANS
	    I=ID_LOW(K)
	    J=ID_UP(K)
	    OMEGA(I,J)=OMEGA_TABLE(ID_INDX(K)+NUM_TVALS-1)
	  END DO
	ELSE
	  L=2
	  DO WHILE(TEMP .GT. T_TABLE(L))
	    L=L+1
	  END DO
	  T1=DLOG(T_TABLE(L)/T_TABLE(L-1))
	  DO K=1,NUM_TRANS
	    I=ID_LOW(K)
	    J=ID_UP(K)
!
! NB: ID_INDX(K) refers to OMEGA for T_TABLE(1). Thus
! NB: ID_INDX(K)+L-1 refers to OMEGA for T_TABLE(L).
!
	    ALPHA=OMEGA_TABLE(ID_INDX(K)+L-1)/
	1                    OMEGA_TABLE(ID_INDX(K)+L-2)
	    ALPHA=DLOG(ALPHA)/T1
	    OMEGA(I,J)=OMEGA_TABLE(ID_INDX(K)+L-2)*
	1                         (TEMP/T_TABLE(L-1))**ALPHA
	    dln_OMEGA_dlnT(I,J)=ALPHA
	  END DO
	END IF
1000	CONTINUE
!
! Use approximate formula (only for bound-bound collisions) for 
! non-tabulated collisions.
!
! For ZION > 1 we use the approximate expression given by Mihalas on
! page 133 (2nd ed).
! For ZION=1 we use the expression given by Auer and Mihals,
! 1973, ApJ, 184, 151
!
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(I,J,FL,X,IFAIL,EX_E1X,GBAR,dln_GBAR_dlnT,G1,G2)
	DO I=1,NLEV
	  DO J=I+1,NLEV
	    IF(OMEGA(I,J) .EQ. 0.0D0 .AND. EIN_A(I,J) .NE. 0.0D0)THEN
	      FL=EDGE(I)-EDGE(J)
	      X=HDKT*FL/TEMP
	      IFAIL=0
	      EX_E1X=EX_E1X_FUN(X,IFAIL)
	      IF(IFAIL .NE. 0)THEN
	        WRITE(LUER,*)'Error in OMEGA_GEN when computing EX_E1X'
	        WRITE(LUER,*)'FILE_NAME=',FILE_NAME
	        WRITE(LUER,*)I,J,X
	        STOP
	      END IF
	      IF(ZION .EQ. 1)THEN
	        IF(X .LE. 14.0D0)THEN
	          GBAR=0.276D0*EX_E1X
	          dln_GBAR_dlnT=(1.0D0/EX_E1X-X)
	        ELSE
	          GBAR=0.066D0*(1.0D0+1.5D0/X)/SQRT(X)
	          dln_GBAR_dlnT=-0.5D0+0.099D0/X/GBAR
	        END IF
	      ELSE
	        G1=0.2D0
	        IF( SAME_N(LEVELNAME(I),LEVELNAME(J)) )G1=0.7D0
	        G2=0.276D0*EX_E1X
	        IF(G1 .GT. G2)THEN
	          GBAR=G1
	          dln_GBAR_dlnT=0.0D0
	        ELSE
	          GBAR=G2
	          dln_GBAR_dlnT=(1.0D0/EX_E1X-X)
	        END IF
	      END IF
	      OMEGA(I,J)=47.972D0*OMEGA_SCALE*GBAR*EIN_A(I,J)*STAT_WT(I)/FL
	      dln_OMEGA_dlnT(I,J)=dln_GBAR_dlnT
	    ELSE IF(OMEGA(I,J) .EQ. 0.0D0)THEN
	      OMEGA(I,J)=OMEGA_SET
	      dln_OMEGA_dlnT(I,J)=0.0D0
	    END IF
	  END DO
	END DO
!
! Compute the collisional ionization rate.
!
! The approximate ionization formulae, using the photionization cross
! section, is from Mihalas (178, p133,134).
!  GBAR=0.1 if Z = 1.
!  GBAR=0.2 if Z = 2.
!  GBAR=0.3 if Z >= 3.
!
! Note that the photionization cross-section is a factor of 10**10 to large.
! This must be compensated for in the constant at the beginning of the
! expression.
!
! Constant is 3.23*0.3/8.63E-08
!
	PHOT_CONST=0.1D0*MIN(ZION,THREE)*3.23D0/8.63D-08
	PHOT_ID=1
	RET_EDGE_CROSS=.TRUE.
	IF(LAST_RD_FILE .NE. FILE_NAME)THEN
	  IF(ALLOCATED(PHOT_CROSS))THEN
	    I=SIZE(PHOT_CROSS)
	    IF(I .LT. NLEV)THEN
	      DEALLOCATE(PHOT_CROSS)
	      ALLOCATE(PHOT_CROSS(NLEV))
	    END IF
	  ELSE
	    ALLOCATE(PHOT_CROSS(NLEV))
	  END IF
	  CALL SUB_PHOT_GEN(ID,PHOT_CROSS,ZERO,EDGE,NLEV,PHOT_ID,RET_EDGE_CROSS)
	END IF
	DO I=1,NLEV
	  IF(OMEGA(I,I) .EQ. 0.0D0)THEN
	    OMEGA(I,I)=PHOT_CONST*PHOT_CROSS(I)*STAT_WT(I)*TEMP/EDGE(I)
	    dln_OMEGA_dlnT(I,I)=1.0D0
	  END IF
	END DO
	LAST_RD_FILE=FILE_NAME
!
	RETURN
	END
