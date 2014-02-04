C
C Routine to increment the photoionization and recombination rates
C for an arbitrary ion. The bound-free cooling rate (in ergs/cm**3/s)
C is also computed. The Free-Free cooling rate is computed under the
C assumption that is is hydrogenic, and the ion has charge ZHYD.
C
C Separate quadrature weights are passed to compute the cooling rate.
C
C
	SUBROUTINE PRRR_SL_V4(PR,RR,BFCR,FF,WSE,WCR,
	1                     HN,HNST,NLEV,ZHYD,
	1                     DI,DIST,N_DI,
	1                     PHOT_ID,ION_LEV,ED,T,
	1                     JREC,JPHOT,JREC_CR,JPHOT_CR,BPHOT_CR,
	1                     NU_CONT,INIT_ARRAYS,ND)
	IMPLICIT NONE
C
C Altered 14-May-2001 : Bug fixed. Arrays were not being initialized
C                         correctly. Using continuum bands, ML may not be
C                         one on first call. Replaced ML by INIT_ARRAYS.
C                         Changed to V4.
C Altered 25-May-1996 : DIM_LIM removed (now use dynamic memory allocation for
C                         GFF_VAL)
C Altered 29-Sep-1995 : DI,DIST inserted to allow treatment of ionizations
C                         to multiple final states without the need of
C                         separate LTE population for each target level in
C                         the final ion.
C                       Call changed extensively. Now version V2
C                       FLAG deleted as testing of whether to initialize
C                         arrays can be done using PHOT_ID.
C
C Created 23-Sep-87 - Based on PRRRCOOLGEN_V3
C
C
C CONSTANTS FOR OPACITY ETC.
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
	INTEGER NLEV			!Number of levels in species
	INTEGER N_DI			!Number of levels in final ion
	INTEGER ND			!Number of depth points.
	INTEGER PHOT_ID		!Photoionization ID
	INTEGER ION_LEV		!Final (destination) level in ion.
C
	REAL*8 PR(NLEV,ND),RR(NLEV,ND),BFCR(NLEV,ND),FF(ND)
	REAL*8 HN(NLEV,ND),HNST(NLEV,ND),WSE(NLEV,ND),WCR(NLEV,ND)
	REAL*8 DI(N_DI,ND),DIST(N_DI,ND)
C
	REAL*8 ED(ND),T(ND)
	REAL*8 JREC(ND)
	REAL*8 JPHOT(ND)
	REAL*8 JREC_CR(ND)
	REAL*8 JPHOT_CR(ND)
	REAL*8 BPHOT_CR(ND)
	REAL*8 NU_CONT
C
	LOGICAL INIT_ARRAYS	        !Used to signify initialization
C
	INTEGER I,J
	REAL*8 T2,A1,TMP_HNST,B_RAT
	REAL*8 H,ZHYD,CHIBF,CHIFF,HDKT,TWOHCSQ
C
C Dynamic memory allocation for free-free gaunt factor as a function of depth.
C
	REAL*8 GFF_VAL(ND)
C
C 4PI*1.0E-10 (R scaling) Note that ordering is important or get underflow.
C FQW is approximately 10^15.
C
	H=6.6261965D-12					!ergs/s (*1.0E+15 due to *nu)
C
C If ML=1 and and PHOT_ID .EQ. 1 then initialize all arrays. This routine
C should be called first for ionizations to the ground state.
C
	IF(INIT_ARRAYS .AND. PHOT_ID .EQ. 1)THEN
	  PR(:,:)=0.0D0				!NLEV,ND
	  RR(:,:)=0.0D0				!NLEV,ND
	  BFCR(:,:)=0.0D0			!NLEV,ND
	  FF(:)=0.0D0				!ND
	END IF
C
C Note that JREC     = Int [ (2hv^3/c^2 +J) exp(-hv/kT)/v dv ]
C           JREC_CR  = Int [ (2hv^3/c^2 +J) exp(-hv/kT)   dv ]
C           JPHOT    = Int [ J/v dv]
C           JPHOT_CR = Int [ J dv]
C
C Since BFCR = Int (nu-edge)/nu, J?_CR is associated with WSE in the expression
C for BFCR.
C
	DO J=1,ND
	  B_RAT=(DI(ION_LEV,J)/DIST(ION_LEV,J))*(DIST(1,J)/DI(1,J))
	  DO I=1,NLEV
	    IF(WSE(I,J) .NE. 0)THEN
	      TMP_HNST=HNST(I,J)*B_RAT
	      PR(I,J)=PR(I,J)+WSE(I,J)*HN(I,J)*JPHOT(J)
	      RR(I,J)=RR(I,J)+WSE(I,J)*TMP_HNST*JREC(J)
	      BFCR(I,J)=BFCR(I,J)+
	1          ( TMP_HNST*(WCR(I,J)*JREC(J)+WSE(I,J)*JREC_CR(J))
	1              -HN(I,J)*(WCR(I,J)*JPHOT(J)+WSE(I,J)*JPHOT_CR(J)) )*H
	    END IF
	  END DO
	END DO
C
C Now compute Free-Free cooling.
C Compute free-free gaunt factors. Replaces call to GFF in following DO loop.
C
	T2=1.256637E-09*ZHYD*ZHYD*CHIFF/(NU_CONT**3)
	CALL GFF_VEC(GFF_VAL,NU_CONT,T,ZHYD,ND)
	DO J=1,ND
	  A1=EXP(-HDKT*NU_CONT/T(J))
	  FF(J) =FF(J)+T2*ED(J)*DI(ION_LEV,J)/SQRT(T(J))*(1.0D0-A1)
	1       *GFF_VAL(J)*( BPHOT_CR(J)-JPHOT_CR(J) )
	END DO
C
	RETURN
	END
