!
! Subroutine to compute the quadrature weights for the statistical
! equilibrium equations. These quadrature weight now have to be multplied
! by FQW/NU before use. This change was made to allow for a fixed continuum
! photioization cross-section.
!
! This routine is specifically designed for the handling of super levels.
! That is, we treat the process in a large atom but assume that the populations
! can be described by a smaller set of levels.
!
! Routine also handles level dissolution.
!
! Notation:
!
!         We use _F to denote populations and variables for the FULL atom,
!            with all terms and levels treated separately.
!	  We use _S to denote populations and variables for the SMALL model
!            atom, with many terms and levels treated as one (i.e using
!            SUPER levels).
!
	SUBROUTINE QUAD_MULTI_V9(WSE_S,dWSE_SdT,WSE_CR_S,
	1                       HNST_S,dlnHNST_S_dlnT,N_S,
	1                       HNST_F_ON_S,EDGE_F,N_F,
	1                       F_TO_S_MAPPING,NU_CONT,T,ND,
	1                       COMPUTE_BA,FIXED_T,LAST_ITERATION,
	1                       DESC,ZION,PHOT_ID,ID)
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
	EXTERNAL SUB_PHOT_GEN
!
! Altered 23-Oct-2016 - Bug fix (collision ionization with dissolution) and introduced
!                          PHOT_DIS_PARAMETER.
! Altered 04-Oct-2016 - Changed to V9
!                       Only compute WSE_CR_S(I_S,J) and dWSE_SdT on LAST iteration, or
!                          when T is variable.
! Altered 05-Apr-2011 - Changed to V8.
!                       HNST_F_ON_S (rather than HNST_F) is passed in call.
!                       HNST_F/HNST_S replaced by HNST_F_ON_S - done to faciliate
!                         modifications allowing lower temperaturs.
!                       Most of editing done early 2011 
! Altered 15-Dec-1997 - MOD_LEV_DIS_BLK replaces include file. Level
!                         dissolution can be switched off completely.
! Altered 05-Sep-1997 - Option to assume that continuum cross-sections have
!                         not altered since the last call. Should result
!                         in CPU time. The definitions of WSE etc have changed
!                         which means that the following routines also need
!                         changing:
!                                   EVALSE
!                         As call changed, now V5.
!                  
! Altered 20-Sep-1996 - Extensive changes to allow SUB_PHOT to be called.
!                       Changes designed to improve speed and vectorization.
!                       Level dissolution effects directy incorporated.
!                       (As extensive changes called _V4, 13-Dec-1996)
! Altered 28-May-1996 - Now call PHOT_GEN_BLEND_V2
! Altered 08-Jun-1995 - WSE_CR_S installed.
! Created 15-May-1995 - Based on QUADGEN_V4
!
	INTEGER ID
	INTEGER N_S,N_F,ND
	REAL*8 WSE_S(N_S,ND)
	REAL*8 dWSE_SdT(N_S,ND)
	REAL*8 WSE_CR_S(N_S,ND)
	REAL*8 HNST_S(N_S,ND)
	REAL*8 dlnHNST_S_dlnT(N_S,ND)
!
	REAL*8 HNST_F_ON_S(N_F,ND)
	REAL*8 EDGE_F(N_F)			!In 10^15 Hz
	INTEGER F_TO_S_MAPPING(N_F)
	REAL*8 T(ND)
!
	REAL*8 NU_CONT
	REAL*8 ZION
	CHARACTER*(*) DESC
	INTEGER PHOT_ID
!
	REAL*8 YDIS(ND)		!Constant for computing level dissolution/
	REAL*8 XDIS(ND)		!Constant for computing level dissolution/
	REAL*8 DIS_CONST(N_F)	!Constant appearing in dissolution formula.
	REAL*8 ALPHA_VEC(N_F)
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 OPLIN,EMLIN
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL COMPUTE_BA,FIXED_T,LAST_ITERATION
!
! Local Variables,
!
	INTEGER I_S,I_F,J
	REAL*8 T1,T2,T3,ZION_CUBED,NEFF,FOUR_PI_D_H
	LOGICAL DO_ALL
!
! NB: WSE_OLD=WSE*FQW/NU
!     dWSEdT_OLD=dWSEdT*FQW/NU
!     WSE_CR_S_OLD=(NU*WSE+WSE_CR_S)*FQW/NU
!
! The factor of DEX(-10) in FOUR_PI_D_H is due to the definition of the
! cross-section in SUB_GEN_PHOT
! which is DEX(10) times the photoionization cross section so that
! CHI*R is constant.
!
! Note FOUR_PI_D_H differs by 10^-15 from original constant in QUADGEN because
! FQW has C units of Hz, not 10^15 Hz.
!
	FOUR_PI_D_H=1.0D0/5.27296D-03                    !1.8965D+02		!4*PI/H*DEX(-10)*DEX(-15)
!
	WSE_S(:,:)=0.0D0
	WSE_CR_S(:,:)=0.0D0
	dWSE_SdT(:,:)=0.0D0
!
! Get photoionization cross-sections for all levels. The first call returns
! the threshold cross-section when NU < EDGE.
!
	IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	  CALL SUB_PHOT_GEN(ID,ALPHA_VEC,NU_CONT,EDGE_F,N_F,PHOT_ID,L_TRUE)
	ELSE
	  CALL SUB_PHOT_GEN(ID,ALPHA_VEC,NU_CONT,EDGE_F,N_F,PHOT_ID,L_FALSE)
	END IF
!
! DIS_CONST is the constant K appearing in the expression for level dissolution.
! A negative value for DIS_CONST implies that the cross-section is zero.
!
	DIS_CONST(1:N_F)=-1.0D0
	IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	  ZION_CUBED=ZION*ZION*ZION
	  DO I_F=1,N_F
	    IF(NU_CONT .LT. EDGE_F(I_F) .AND. ALPHA_VEC(I_F) .NE. 0 .AND. NU_CONT .GT. 0.8D0*EDGE_F(I_F))THEN
	      NEFF=SQRT(3.289395D0*ZION*ZION/(EDGE_F(I_F)-NU_CONT))
	      IF(NEFF .GT. 2*ZION)THEN
	        T1=MIN(1.0D0,16.0D0*NEFF/(1+NEFF)/(1+NEFF)/3.0D0)
	         DIS_CONST(I_F)=( T1*ZION_CUBED/(NEFF**4) )**1.5D0
	      END IF
	    END IF
	  END DO
	END IF
!
! Compute dissolution vectors that are independent of level.
!
	IF(MOD_DO_LEV_DIS)THEN
	  DO J=1,ND
	    YDIS(J)=1.091D0*(X_LEV_DIS(J)+4.0D0*(ZION-1)*A_LEV_DIS(J))*
	1               B_LEV_DIS(J)*B_LEV_DIS(J)
	    XDIS(J)=B_LEV_DIS(J)*X_LEV_DIS(J)
	  END DO
	END IF
!
	DO_ALL=.FALSE.
	IF(COMPUTE_BA)DO_ALL=.TRUE.
	IF(FIXED_T)DO_ALL=.FALSE.
	IF(LAST_ITERATION)DO_ALL=.TRUE.
!
! We have to loop over depth (rather than frequency) because of the
! FULL to SUPER level mapping.
!
	IF(DO_ALL)THEN
	  DO I_F=1,N_F
	    I_S=F_TO_S_MAPPING(I_F)
	    IF(NU_CONT .GE. EDGE_F(I_F))THEN
	      T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	      DO J=1,ND
	        WSE_S(I_S,J)=WSE_S(I_S,J) + T1*HNST_F_ON_S(I_F,J)
	        WSE_CR_S(I_S,J)=WSE_CR_S(I_S,J) - EDGE_F(I_F)*T1*HNST_F_ON_S(I_F,J)
	        dWSE_SdT(I_S,J)=dWSE_SdT(I_S,J) - T1*HNST_F_ON_S(I_F,J)*
	1          (dlnHNST_S_dlnT(I_S,J)+1.5D0+HDKT*EDGE_F(I_F)/T(J))/T(J)
	      END DO
!
! We only allow for level dissolutions when the ionizations are occurring to
! the ground state.
!
	    ELSE IF(DIS_CONST(I_F) .GE. 0.0D0)THEN
	      T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	      DO J=1,ND
	        T2=7.782D0+XDIS(J)*DIS_CONST(I_F)
	        T3=T2/(T2+YDIS(J)*DIS_CONST(I_F)*DIS_CONST(I_F))
	        IF(T2 .GT. PHOT_DIS_PARAMETER)THEN
	          T3=T1*T3
	          WSE_S(I_S,J)=WSE_S(I_S,J) + T3*HNST_F_ON_S(I_F,J)
	          WSE_CR_S(I_S,J)=WSE_CR_S(I_S,J) - EDGE_F(I_F)*T3*HNST_F_ON_S(I_F,J)
	          dWSE_SdT(I_S,J)=dWSE_SdT(I_S,J) - T3*HNST_F_ON_S(I_F,J)*
	1            (dlnHNST_S_dlnT(I_S,J)+1.5D0+HDKT*EDGE_F(I_F)/T(J))/T(J)
	        END IF
	      END DO
	    END IF
	  END DO
	ELSE
	  DO I_F=1,N_F
	    I_S=F_TO_S_MAPPING(I_F)
	    IF(NU_CONT .GE. EDGE_F(I_F))THEN
	      T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	      DO J=1,ND
	        WSE_S(I_S,J)=WSE_S(I_S,J) + T1*HNST_F_ON_S(I_F,J)
	      END DO
!
! We only allow for level dissolutions when the ionizations are occurring to
! the ground state.
!
	    ELSE IF(DIS_CONST(I_F) .GE. 0.0D0)THEN
	      T1=FOUR_PI_D_H*ALPHA_VEC(I_F)
	      DO J=1,ND
	        T2=7.782D0+XDIS(J)*DIS_CONST(I_F)
	        T3=T2/(T2+YDIS(J)*DIS_CONST(I_F)*DIS_CONST(I_F))
	        IF(T2 .GT. PHOT_DIS_PARAMETER)THEN
	          T3=T1*T3
	          WSE_S(I_S,J)=WSE_S(I_S,J) + T3*HNST_F_ON_S(I_F,J)
	        END IF
	      END DO
	    END IF
	  END DO
	END IF
!  
	RETURN
	END
