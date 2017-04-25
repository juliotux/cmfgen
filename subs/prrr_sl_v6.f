!
! Routine to increment the photoionization and recombination rates
! for an arbitrary ion. The bound-free cooling rate (in ergs/cm**3/s)
! is also computed. The Free-Free cooling rate is computed under the
! assumption that is is hydrogenic, and the ion has charge ZHYD.
!
! Separate quadrature weights are passed to compute the cooling rate.
!
!
	SUBROUTINE PRRR_SL_V6(PR,RR,BFCR,FF,WSE,WCR,
	1                     HN,HNST,NLEV,ZHYD,
	1                     DI,LOG_DIST,N_DI,
	1                     PHOT_ID,ION_LEV,ED,T,
	1                     JREC,JPHOT,JREC_CR,JPHOT_CR,BPHOT_CR,
	1                     NU,NU_CONT,INIT_ARRAYS,ND)
	IMPLICIT NONE
!
! Altered 23-Jun-2015 - Added H- free-free cooling.
! Altered 20-Oct-2011 - Now sum up all in ion levels for FF. Only do this when PHOT_ID=1
! Altered 05-Apr-2011 - Changed to V6.
!                       LOG_DIST (rather than dwHNST_F) is passed in call.
!                         Modifications done to allow lower temperaturs.
!                         Primary editing done 25-Jan-2011
! Altered 25-Jan-2010 : Bug fixed with previous alteration.
! Created 29-Nov-2010 : Based on PRRR_SL_V5
!                         LOG_DIST installed to prevent crashing caused by high ionization
!                         stages at low temperatures.
! Altered 03-Mar-2004 : NU installed in CALL. Computation of FF cooling revised.
! Altered 14-May-2001 : Bug fixed. Arrays were not being initialized
!                         correctly. Using continuum bands, ML may not be
!                         one on first call. Replaced ML by INIT_ARRAYS.
!                         Changed to V4.
! Altered 25-May-1996 : DIM_LIM removed (now use dynamic memory allocation for
!                         GFF_VAL)
! Altered 29-Sep-1995 : DI,DIST inserted to allow treatment of ionizations
!                         to multiple final states without the need of
!                         separate LTE population for each target level in
!                         the final ion.
!                       Call changed extensively. Now version V2
!                       FLAG deleted as testing of whether to initialize
!                         arrays can be done using PHOT_ID.
!
! Created 23-Sep-87 - Based on PRRRCOOLGEN_V3
!
!
! CONSTANTS FOR OPACITY ETC.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER NLEV			!Number of levels in species
	INTEGER N_DI			!Number of levels in final ion
	INTEGER ND			!Number of depth points.
	INTEGER PHOT_ID		!Photoionization ID
	INTEGER ION_LEV		!Final (destination) level in ion.
!
	REAL*8 PR(NLEV,ND),RR(NLEV,ND),BFCR(NLEV,ND),FF(ND)
	REAL*8 HN(NLEV,ND),HNST(NLEV,ND),WSE(NLEV,ND),WCR(NLEV,ND)
	REAL*8 DI(N_DI,ND),LOG_DIST(N_DI,ND)
!
	REAL*8 ED(ND),T(ND)
	REAL*8 JREC(ND)
	REAL*8 JPHOT(ND)
	REAL*8 JREC_CR(ND)
	REAL*8 JPHOT_CR(ND)
	REAL*8 BPHOT_CR(ND)
	REAL*8 NU
	REAL*8 NU_CONT
!
	LOGICAL INIT_ARRAYS	        !Used to signify initialization
!
	INTEGER I,J
	REAL*8 POP_SUM,T2,A1,TMP_HNST
	REAL*8 JB_RAT, JC_RAT
	REAL*8 H,ZHYD,CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Dynamic memory allocation for free-free gaunt factor as a function of depth.
!
	REAL*8 GFF_VAL(ND)
!
! 4PI*1.0E-10 (R scaling) Note that ordering is important or get underflow.
! FQW is approximately 10^15.
!
	H=6.6261965D-12					!ergs/s (*1.0E+15 due to *nu)
!
! If ML=1 and and PHOT_ID .EQ. 1 then initialize all arrays. This routine
! should be called first for ionizations to the ground state.
!
	IF(INIT_ARRAYS .AND. PHOT_ID .EQ. 1)THEN
	  PR(:,:)=0.0D0				!NLEV,ND
	  RR(:,:)=0.0D0				!NLEV,ND
	  BFCR(:,:)=0.0D0			!NLEV,ND
	  FF(:)=0.0D0				!ND
	END IF
!
! Note that JREC     = Int [ (2hv^3/c^2 +J) exp(-hv/kT)/v dv ]
!           JREC_CR  = Int [ (2hv^3/c^2 +J) exp(-hv/kT)   dv ]
!           JPHOT    = Int [ J/v dv]
!           JPHOT_CR = Int [ J dv]
!
! Since BFCR = Int (nu-edge)/nu, J?_CR is associated with WSE in the expression
! for BFCR.
!
	DO J=1,ND
	  IF(JREC(J) .GT. 0.0D0)THEN
	    JB_RAT=LOG(DI(ION_LEV,J)/DI(1,J))-LOG_DIST(ION_LEV,J)+LOG_DIST(1,J)
	    JB_RAT=EXP(LOG(JREC(J))+JB_RAT)
	  ELSE
	    JB_RAT=0.0D0
	  END IF
	  IF(JREC_CR(J) .GT. 0.0D0)THEN
	    JC_RAT=LOG(DI(ION_LEV,J)/DI(1,J))-LOG_DIST(ION_LEV,J)+LOG_DIST(1,J)
	    JC_RAT=EXP(LOG(JREC_CR(J))+JC_RAT)
	  ELSE
	    JC_RAT=0.0D0
	  END IF
	  DO I=1,NLEV
	    IF(WSE(I,J) .NE. 0)THEN
	      PR(I,J)=PR(I,J)+WSE(I,J)*HN(I,J)*JPHOT(J)
	      RR(I,J)=RR(I,J)+WSE(I,J)*HNST(I,J)*JB_RAT
	      BFCR(I,J)=BFCR(I,J)+
	1          ( HNST(I,J)*(WCR(I,J)*JB_RAT+WSE(I,J)*JC_RAT)
	1              -HN(I,J)*(WCR(I,J)*JPHOT(J)+WSE(I,J)*JPHOT_CR(J)) )*H
	    END IF
	  END DO
	END DO
!
! Compute Free-Free cooling.
!
	IF(ZHYD .EQ. 0.0D0)THEN
!
! This is for H-. We use GFF_VAL as a temporary storage for the ground state
! population of neutral hydrogen.
!
	  GFF_VAL(1:ND)=DI(1,1:ND)
	  CALL DO_HMI_FF_COOL(FF,GFF_VAL,ED,T,BPHOT_CR,JPHOT_CR,NU_CONT,ND)
	ELSE IF(ION_LEV .EQ. 1)THEN
!
! Compute free-free gaunt factors. Replaces call to GFF in following DO loop.
!
	  CALL GFF_VEC(GFF_VAL,NU_CONT,T,ZHYD,ND)
!
! The opacity is evaluated at NU_CONT. However, the stimulated emission term
! should be evaluated at NU, since we correct CHI and ETA for the change
! in freqency but assuming a constant cross-section.
!
! The constant in T2 is 4PI x 1.0E-10.
!
	  T2=1.256637061D-09*ZHYD*ZHYD*CHIFF/(NU_CONT**3)
	  DO J=1,ND
	    POP_SUM=SUM(DI(:,J))
	    A1=EXP(-HDKT*NU/T(J))
	    FF(J) =FF(J)+T2*ED(J)*POP_SUM/SQRT(T(J))*(1.0D0-A1)
	1           *GFF_VAL(J)*( BPHOT_CR(J)-JPHOT_CR(J) )
	  END DO
	END IF
!
	RETURN
	END
