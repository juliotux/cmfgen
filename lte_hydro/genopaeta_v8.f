C
C Subroutine to compute the contribution to the opacity AND emissivity
C by FREE-FREE and BOUND-FREE processes for a general ion. The
C contribution is added directly to the opacity CHI and emissivity
C ETA.
C
	SUBROUTINE GENOPAETA_V8(ID,CHI,ETA,NU,
	1              HN,HNST,EDGE,GION,ZION,N,
	1              DI,DIST,N_DI,PHOT_ID,ION_LEV,
	1              ED,T,EMHNUKT,
	1              IONFF,ND,LST_DEPTH_ONLY)
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
C
C Altered 27-May-2005 - Added check to get minimum frequency, in case levels
C                         are out of order.
C Altered 22-Feb-2000 - ID added to call.
C                       SUB_PHOT_GEN deleted.
C Altered 15-Dec-1997 - MOD_LEV_DIS_BLK replaces include file. Level
C                         dissolution can be switched off completely.
C Altered 17-Sep-1996 - Extensive modifications made to make the routine
C                         more vectorizeable. PHOT_GEN_BLEND inlined
C                         Computation of all photioization cross-sections
C                          removed to a separate subroutine.
C                       As major changes called _V6 (13-Dec-1996)
C
C Altered 26-Aug-1996 - COR_FAC included to improve speed.
C Altered 24-May-1996 - DIM_LIM and ND_MAX removed using F90.
C Altered 28-Sep-1995 - Extensive modifications. Call altered.
C
C Altered 07-Jul-1995 - Bug fix. CROSS(ND) was being incorrectly computed
C                       when LST_DEPTH_ONLY was set to TRUE.
C                       PGOT_GEN_BLEND call changed to _V2.
C
C Altered 08-Jun-1995 - PHOT_GEN_BLEND is now only called when NU < EDGE.
C Altered 02-Jun-1995 - LST_DEPTH_ONLY variable installed. Useful to save
C                         compuatational time when dTDdR is being computed,
C                         and hence when we only require opacities/emissivities
C                         at inner boundary. (Change from _V3 to _V4).
C Altered 12-Aug-1994 - Contribution to opacities by continous level
C                       dissolution included. We do not consider the
C                       B ionizations, as this requires special treatment.
C Altered 07-Jan-1991 - Call GFF_VEC relaces GFF function.
C Altered 30-Oct-1989 - EXTERNAL specification installed for UNIX
C                       compatibility.
C Created 20-Mar-1989 - Based on OPAGEN and CHIGEN
C
	INTEGER ID,N,N_DI,ND
	LOGICAL IONFF,LST_DEPTH_ONLY
C
C Constants for opacity etc.
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
	REAL*8 CHI(ND)			!Opacity
	REAL*8 ETA(ND)			!Emissivity
C
C Large Model Atom Populations.
C
	REAL*8 HN(N,ND),HNST(N,ND)
	REAL*8 EDGE(N)
C
C Ion populations. These populations should refer to the small model atoms.
C (i.e. the model atom with super levels_
C
	REAL*8 DI(N_DI,ND),DIST(N_DI,ND)
C
	REAL*8 T(ND)			!Temperature (K)
	REAL*8 ED(ND)			!Electron density
	REAL*8 EMHNUKT(ND)		!EXP(-hv/kT)
	REAL*8 NU			!Frequency (10^15 Hz)
	REAL*8 ZION			!Charge on resultiong ion.
	REAL*8 GION			!Charge on resultiong ion.
C
	INTEGER PHOT_ID		!Photoionization ID (path)
	INTEGER ION_LEV		!Target level for ionizations in ion.
C
C Vectors to save computational effort.
C
	REAL*8 GFF_VAL(ND)		!g(ff) as a function of depth
	REAL*8 COR_FAC(ND)		!Factor to convert HNST for ION_LEV
	REAL*8 YDIS(ND)			!Constant for computing level dissolution/
	REAL*8 XDIS(ND)			!Constant for computing level dissolution/
	REAL*8 DIS_CONST(N)		!Constant appearing in dissolution formula.
	REAL*8 ALPHA_VEC(N)		!Photionization cross-section
	REAL*8 TMP_CHI(N)		!Photionization cross-section
	REAL*8 TMP_ETA(N)		!Photionization cross-section
C
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
C
C Local constants.
C
	INTEGER I,K,K_ST,ND_LOC,NO_NON_ZERO_PHOT
	REAL*8 ALPHA,TCHI1,TETA1,TETA2
	REAL*8 T1,T2,ZION_CUBED,NEFF
C
	CALL TUNE(1,'GENOPA')
C
C ND_LOC indicates the number of depth points we are going to compute the
C opacity at. K_ST indicates the depth point to start, and is either 1, or ND.
C
	IF(LST_DEPTH_ONLY)THEN
	  ND_LOC=1
	  K_ST=ND
	ELSE
	  ND_LOC=ND
	  K_ST=1
	END IF
C
C Add in free-free contribution. We assume that only ground state free-free
C is important. We will generally include the first excited anyway,
C when we compute opacities contributed by ionizations to excited levevels.
C
	IF(IONFF)THEN
C
C Compute free-free gaunt factors. Replaces call to GFF in following DO loop.
C
	  IF(LST_DEPTH_ONLY)THEN
	    CALL GFF_VEC(GFF_VAL(ND),NU,T(ND),ZION,ND_LOC)
	    IF(ION_LEV .EQ. 1)THEN
	      CALL FF_RES_GAUNT(GFF_VAL(ND),NU,T(ND),ID,GION,ZION,ND_LOC)
	    END IF
	  ELSE
	    CALL GFF_VEC(GFF_VAL,NU,T,ZION,ND)
	    IF(ION_LEV .EQ. 1)THEN
	      CALL FF_RES_GAUNT(GFF_VAL,NU,T,ID,GION,ZION,ND)
	    END IF
	  END IF
C
	  TCHI1=CHIFF*ZION*ZION/(NU*NU*NU)
	  TETA1=CHIFF*ZION*ZION*TWOHCSQ
	  DO K=K_ST,ND
	    ALPHA=ED(K)*DI(ION_LEV,K)*GFF_VAL(K)/SQRT(T(K))
	    CHI(K)=CHI(K)+TCHI1*ALPHA*(1.0D0-EMHNUKT(K))
	    ETA(K)=ETA(K)+TETA1*ALPHA*EMHNUKT(K)
	  END DO
	END IF
C
C 
C Now add in BOUND-FREE contributions. We first compute vectors which can
C decrease the compuation time.
C
C NB: A clearer way of writing the expressions for ETA and CHI is
C where TMP_HNST is the LTE population defined by the actual population of
C the destination (target) level.
C
C TMP_HNST=HNST(I,K)*(DI(ION_LEV,K)/DIST(ION_LEV,K))*(DIST(1,K)/DI(1,K))
C CHI(K)=CHI(K)+ALPHA*(HN(I,K)-TMP_HNST*EMHNUKT(K))
C ETA(K)=ETA(K)+TETA2*TMP_HNST*EMHNUKT(K)
C
C In case some levels are out of order, we get make sure we get the minimum edge
C frequency.
C
	T1=MINVAL(EDGE(1:N))
	IF(NU .GE. T1)THEN
	  DO K=K_ST,ND
	    IF(DIST(ION_LEV,K) .NE. 0.0D0 .AND. DI(1,K) .NE. 0)THEN
	      COR_FAC(K)=(DI(ION_LEV,K)/DIST(ION_LEV,K))*
	1                   (DIST(1,K)/DI(1,K))*EMHNUKT(K)
	    ELSE
	       COR_FAC=0.0D0
	    END IF
	  END DO
	END IF
C
C Compute the photo-ionization cross-sections for all levels.
C
	IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	  CALL SUB_PHOT_GEN(ID,ALPHA_VEC,NU,EDGE,N,PHOT_ID,L_TRUE)
	ELSE
	  CALL SUB_PHOT_GEN(ID,ALPHA_VEC,NU,EDGE,N,PHOT_ID,L_FALSE)
	END IF
	NO_NON_ZERO_PHOT=COUNT(ALPHA_VEC .GT. 0.0D0)
	IF(NO_NON_ZERO_PHOT .EQ. 0)RETURN
C
C DIS_CONST is the constant K appearing in the expression for level dissolution.
C A negative value for DIS_CONST implies that the cross-section is zero.
C
	DIS_CONST(1:N)=-1.0D0
	IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	  ZION_CUBED=ZION*ZION*ZION
	  DO I=1,N
	    IF(NU .LT. EDGE(I) .AND. ALPHA_VEC(I) .NE. 0)THEN
	      NEFF=SQRT(3.289395*ZION*ZION/(EDGE(I)-NU))
	      IF(NEFF .GT. 2*ZION)THEN
	        T1=MIN(1.0D0,16.0D0*NEFF/(1+NEFF)/(1+NEFF)/3.0D0)
	        DIS_CONST(I)=( T1*ZION_CUBED/(NEFF**4) )**1.5D0
	      END IF
	    END IF
	  END DO
	END IF
C
C Compute dissolution vectors that are independent of level.
C
	IF(MOD_DO_LEV_DIS)THEN
	  DO K=K_ST,ND
	    YDIS(K)=1.091*(X_LEV_DIS(K)+4.0D0*(ZION-1)*A_LEV_DIS(K))*
	1                 B_LEV_DIS(K)*B_LEV_DIS(K)
	    XDIS(K)=B_LEV_DIS(K)*X_LEV_DIS(K)
	  END DO
	END IF
C
C 
C Now do the actual Bound-Free computation.
C
C Evaluate the bound-free contributions. If N is small, the inner loop
C is over ND, otherwise the inner loop is over N. For some species
C (e.g. FeIV) N can be large (e.g. 300). The optimal switching point is
C unclear because some photo-ionization cross-sections can be zero,
C and because of the different overheads. It is proably machine dependent.
C
	TETA1=TWOHCSQ*(NU**3)
	IF( NO_NON_ZERO_PHOT .LT. 2*(ND-K_ST+1) )THEN
	  DO I=1,N
	    IF(NU .GE. EDGE(I) .AND. ALPHA_VEC(I) .GT. 0.0D0)THEN
	      TETA2=TETA1*ALPHA_VEC(I)
	      DO K=K_ST,ND
	        CHI(K)=CHI(K)+ALPHA_VEC(I)*(HN(I,K)-HNST(I,K)*COR_FAC(K))
	        ETA(K)=ETA(K)+TETA2*HNST(I,K)*COR_FAC(K)
	      END DO
	    ELSE IF(DIS_CONST(I) .GE. 0)THEN
C
C Add in BOUND-FREE contributions due to level dissolution.
C
	      TETA2=TETA1*ALPHA_VEC(I)
	      DO K=K_ST,ND
	        T1=7.782+XDIS(K)*DIS_CONST(I)
	        T2=T1/(T1+YDIS(K)*DIS_CONST(I)*DIS_CONST(I))
	        CHI(K)=CHI(K)+ALPHA_VEC(I)*T2*(HN(I,K)-HNST(I,K)*EMHNUKT(K))
	        ETA(K)=ETA(K)+TETA2*T2*HNST(I,K)*EMHNUKT(K)
	      END DO
	    END IF		!NU > EDGE
	  END DO		!Depth
C
	ELSE
C
C Inner loop is over level.
C
	  DO K=K_ST,ND
	    TMP_CHI(1:N)=0.0D0
	    TMP_ETA(1:N)=0.0D0
	    DO I=1,N
	      IF(NU .GE. EDGE(I) .AND. ALPHA_VEC(I) .GT. 0.0D0)THEN
	        TMP_CHI(I)=ALPHA_VEC(I)*(HN(I,K)-HNST(I,K)*COR_FAC(K))
	        TMP_ETA(I)=ALPHA_VEC(I)*HNST(I,K)*COR_FAC(K)
	      ELSE IF(DIS_CONST(I) .GE. 0)THEN
C
C Add in BOUND-FREE contributions due to level dissolution.
C
	        T1=7.782+XDIS(K)*DIS_CONST(I)
	        T2=T1/(T1+YDIS(K)*DIS_CONST(I)*DIS_CONST(I))
	        TMP_CHI(I)=ALPHA_VEC(I)*T2*(HN(I,K)-HNST(I,K)*EMHNUKT(K))
	        TMP_ETA(I)=ALPHA_VEC(I)*T2*HNST(I,K)*EMHNUKT(K)
	      END IF
	    END DO		!Level
	    CHI(K)=CHI(K)+SUM(TMP_CHI)
	    ETA(K)=ETA(K)+TETA1*SUM(TMP_ETA)
	  END DO		!Depth
	END IF			!Which inner loop.
C
	CALL TUNE(2,'GENOPA')
	RETURN
	END
