C
C Subroutine to compute the quadrature weights for the statistical
C equilibrium equations. These quadrature weight now have to be multplied
C by FQW/NU before use. This change was made to allow for a fixed continuum
C photioization cross-section.
C
C This routine is specifically designed for the handling of super levels.
C That is, we treat the process in a large atom but assume that the populations
C can be described by a smaller set of levels.
C
C Routine also handles level dissolution.
C
C Notation:
C
C         We use _F to denote populations and variables for the FULL atom,
C            with all terms and levels treated separately.
C	  We use _S to denote populations and variables for the SMALL model
C            atom, with many terms and levels treated as one (i.e using
C            SUPER levels).
C
	SUBROUTINE QUAD_MULTI_V7(WSE_S,dWSE_SdT,WSE_CR_S,
	1                       HNST_S,dlnHNST_S_dlnT,N_S,
	1                       HNST_F,EDGE_F,N_F,
	1                       F_TO_S_MAPPING,NU_CONT,T,ND,
	1                       DESC,ZION,PHOT_ID,ID)
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
	EXTERNAL SUB_PHOT_GEN
C
C Altered 15-Dec-1997 - MOD_LEV_DIS_BLK replaces include file. Level
C                         dissolution can be switched off completely.
C Altered 05-Sep-1997 - Option to assume that continuum cross-sections have
C                         not altered since the last call. Should result
C                         in CPU time. The definitions of WSE etc have changed
C                         which means that the following routines also need
C                         changing:
C                                   EVALSE
C                         As call changed, now V5.
C                  
C Altered 20-Sep-1996 - Extensive changes to allow SUB_PHOT to be called.
C                       Changes designed to improve speed and vectorization.
C                       Level dissolution effects directy incorporated.
C                       (As extensive changes called _V4, 13-Dec-1996)
C Altered 28-May-1996 - Now call PHOT_GEN_BLEND_V2
C Altered 08-Jun-1995 - WSE_CR_S installed.
C Created 15-May-1995 - Based on QUADGEN_V4
C
	INTEGER ID
	INTEGER N_S,N_F,ND
	REAL*8 WSE_S(N_S,ND)
	REAL*8 dWSE_SdT(N_S,ND)
	REAL*8 WSE_CR_S(N_S,ND)
	REAL*8 HNST_S(N_S,ND)
	REAL*8 dlnHNST_S_dlnT(N_S,ND)
C
	REAL*8 HNST_F(N_F,ND)
	REAL*8 EDGE_F(N_F)			!In 10^15 Hz
	INTEGER F_TO_S_MAPPING(N_F)
	REAL*8 T(ND)
C
	REAL*8 NU_CONT
	REAL*8 ZION
	CHARACTER*(*) DESC
	INTEGER PHOT_ID
C
	REAL*8 YDIS(ND)		!Constant for computing level dissolution/
	REAL*8 XDIS(ND)		!Constant for computing level dissolution/
	REAL*8 DIS_CONST(N_F)	!Constant appearing in dissolution formula.
	REAL*8 ALPHA_VEC(N_F)
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 OPLIN,EMLIN
C
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
C
C Local Variables,
C
	INTEGER I_S,I_F,J
	REAL*8 T1,T2,T3,ZION_CUBED,NEFF,FOUR_PI_D_H
C
C NB: WSE_OLD=WSE*FQW/NU
C     dWSEdT_OLD=dWSEdT*FQW/NU
C     WSE_CR_S_OLD=(NU*WSE+WSE_CR_S)*FQW/NU
C
C The factor of DEX(-10) in FOUR_PI_D_H is due to the definition of the
C cross-section in SUB_GEN_PHOT
C which is DEX(10) times the photoionization cross section so that
C CHI*R is constant.
C
C Note FOUR_PI_D_H differs by 10^-15 from original constant in QUADGEN because
C FQW has C units of Hz, not 10^15 Hz.
C
	FOUR_PI_D_H=1.0D0/5.27296D-03                    !1.8965D+02		!4*PI/H*DEX(-10)*DEX(-15)
C
	WSE_S(:,:)=0.0D0
	WSE_CR_S(:,:)=0.0D0
	dWSE_SdT(:,:)=0.0D0
C
C Get photoionization cross-sections for all levels. The first call returns
C the threshold cross-section when NU < EDGE.
C
	IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	  CALL SUB_PHOT_GEN(ID,ALPHA_VEC,NU_CONT,EDGE_F,N_F,PHOT_ID,L_TRUE)
	ELSE
	  CALL SUB_PHOT_GEN(ID,ALPHA_VEC,NU_CONT,EDGE_F,N_F,PHOT_ID,L_FALSE)
	END IF
C
C DIS_CONST is the constant K appearing in the expression for level dissolution.
C A negative value for DIS_CONST implies that the cross-section is zero.
C
	DIS_CONST(1:N_F)=-1.0D0
	IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	  ZION_CUBED=ZION*ZION*ZION
	  DO I_F=1,N_F
	    IF(NU_CONT .LT. EDGE_F(I_F) .AND. ALPHA_VEC(I_F) .NE. 0)THEN
	      NEFF=SQRT(3.289395*ZION*ZION/(EDGE_F(I_F)-NU_CONT))
	      IF(NEFF .GT. 2*ZION)THEN
	        T1=MIN(1.0D0,16.0D0*NEFF/(1+NEFF)/(1+NEFF)/3.0D0)
	         DIS_CONST(I_F)=( T1*ZION_CUBED/(NEFF**4) )**1.5D0
	      END IF
	    END IF
	  END DO
	END IF
C
C Compute dissolution vectors that are independent of level.
C
	IF(MOD_DO_LEV_DIS)THEN
	  DO J=1,ND
	    YDIS(J)=1.091*(X_LEV_DIS(J)+4.0D0*(ZION-1)*A_LEV_DIS(J))*
	1               B_LEV_DIS(J)*B_LEV_DIS(J)
	    XDIS(J)=B_LEV_DIS(J)*X_LEV_DIS(J)
	  END DO
	END IF
C
C We have to loop over depth (rather than frequency) because of the
C FULL to SUPER level mapping.
C
	DO I_F=1,N_F
	  I_S=F_TO_S_MAPPING(I_F)
	  IF(NU_CONT .GE. EDGE_F(I_F))THEN
	    T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	    DO J=1,ND
	      WSE_S(I_S,J)=WSE_S(I_S,J) +
	1        T1*(HNST_F(I_F,J)/HNST_S(I_S,J))
	      WSE_CR_S(I_S,J)=WSE_CR_S(I_S,J) - EDGE_F(I_F)*
	1           T1*(HNST_F(I_F,J)/HNST_S(I_S,J))
	      dWSE_SdT(I_S,J)=dWSE_SdT(I_S,J) -
	1        T1*(HNST_F(I_F,J)/HNST_S(I_S,J))*
	1        (dlnHNST_S_dlnT(I_S,J)+1.5D0+HDKT*EDGE_F(I_F)/T(J))/T(J)
	    END DO
C
C We only allow for level dissolutions when the ionizations are occurring to
C the ground state.
C
	  ELSE IF(DIS_CONST(I_F) .GE. 0)THEN
	    T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	    DO J=1,ND
	      T2=7.782+XDIS(J)*DIS_CONST(I_F)
	      T3=T1*T2/(T2+YDIS(J)*DIS_CONST(I_F)*DIS_CONST(I_F))
	      WSE_S(I_S,J)=WSE_S(I_S,J) +
	1           T3*(HNST_F(I_F,J)/HNST_S(I_S,J))
	      WSE_CR_S(I_S,J)=WSE_CR_S(I_S,J) - EDGE_F(I_F)*
	1           T3*(HNST_F(I_F,J)/HNST_S(I_S,J))
	      dWSE_SdT(I_S,J)=dWSE_SdT(I_S,J) -
	1           T3*(HNST_F(I_F,J)/HNST_S(I_S,J))*
	1          (dlnHNST_S_dlnT(I_S,J)+1.5D0+HDKT*EDGE_F(I_F)/T(J))/T(J)
	    END DO
	  END IF
	END DO
	    
	RETURN
	END
