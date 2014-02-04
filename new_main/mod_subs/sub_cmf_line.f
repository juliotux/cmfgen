!
! Subroutine to compute JBAR and update the variation matrix (BA) for
! lines treated using the solution of the transfer equation in the comoving
! frame without overlap. Routine is based on a section of code in CMFGEN_SUB.F and
! still uses LINEGEN.INC.
!
! Routine was tested and appears to work. This option is becoming obsolete in CMFGEN.
!
	SUBROUTINE SUB_CMF_LINE(SECTION,POPS,CHIL,ETAL,NEG_OPACITY,
	1                    FL,CONT_FREQ,AMASS,EDDINGTON,IMPURITY_CODE,
	1                    EW,CONT_INT,COMPUTE_EW,
	1                    COMPUTE_JEW,LU_JEW,ACCESS_JEW,
	1                    NL,NUP,NT,ND,NC,NP,
	1                    NDEXT,NCEXT,NPEXT,
	1                    NLF,NNM,
	1                    DIAG_INDX,NUM_BNDS)
!
	USE ANG_QW_MOD
	USE CONTROL_VARIABLE_MOD
	USE CMF_SOB_MOD 
	USE LINE_MOD
	USE LINE_VEC_MOD
	USE MOD_CMFGEN
	USE OPAC_MOD
	USE RADIATION_MOD
	USE STEQ_DATA_MOD
 	USE VAR_RAD_MOD
	IMPLICIT NONE
!
! Cleaned : 17-Dec-2004
!
	REAL*8 FL
	REAL*8 CONT_FREQ
	REAL*8 AMASS
	REAL*8 EW
	REAL*8 CONT_INT
	LOGICAL COMPUTE_JEW
	LOGICAL COMPUTE_EW
	INTEGER LU_JEW
	INTEGER ACCESS_JEW
!
	INTEGER NL
	INTEGER NUP
	INTEGER NT
	INTEGER ND,NC,NP
	INTEGER NDEXT,NCEXT,NPEXT
	INTEGER NLF
	INTEGER NNM
	INTEGER DIAG_INDX
	INTEGER NUM_BNDS
	LOGICAL EDDINGTON
	LOGICAL IMPURITY_CODE
	CHARACTER(LEN=*) SECTION
!
	REAL*8 POPS(NT,ND)
	REAL*8 CHIL(ND)
	REAL*8 ETAL(ND)
	LOGICAL NEG_OPACITY(ND)
!
! Constants for opacity etc.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 OPLIN,EMLIN
!
! Local variables.
!
	REAL*8 VB(ND)
	REAL*8 VC(ND)
	REAL*8 T1,T2,T3,T4
	INTEGER I,J,K,L
	INTEGER MNL,MNL_F
	INTEGER MNUP,MNUP_F
	INTEGER MNT
	INTEGER JJ
	INTEGER ID
	INTEGER PHOT_ID
	INTEGER FREQ_INDX
	INTEGER LUER,ERROR_LU
	INTEGER COUNT
	EXTERNAL ERROR_LU
!
	INTEGER GET_DIAG
	INTEGER BNDST
	INTEGER BNDEND
	INTEGER BND_TO_FULL
	LOGICAL LST_DEPTH_ONLY 
	LOGICAL PRINT_DIAGNOSTICS
	LOGICAL COMPUTE_LAM
	LOGICAL FULL_ES
	DATA COUNT/0/
!
! This function takes a band-index and converts it the equivalent index
! in the full matrix. L=BND_TO_FULL(J,K) is equivalent to the statements:
!     IF(NUM_BNDS .EQ. ND)THEN L=J ELSE L=K+J-DIAG_INDX END IF
! The second indice is the equation depth.
!
	BND_TO_FULL(J,K)=(NUM_BNDS/ND)*(DIAG_INDX-K)+K+J-DIAG_INDX
!
! This function computes the index L on BA( , ,?,K) corresponding
! to the local depth variable (i.e that at K). It is equivalent
! to IF (NUM_BNDS .EQ. ND)THEN L=K ELSE L=DIAG END IF
!
	GET_DIAG(K)=(NUM_BNDS/ND)*(K-DIAG_INDX)+DIAG_INDX
!
! These two functions compute the start and end indices when updating
! VJ. eg. we do not wish to update VJ( ,1, ) if we are using the banded
! matrix since this refers to a variable beyond the outer atmosphere.
!
	BNDST(K)=MAX( (NUM_BNDS/ND)*(K-DIAG_INDX)+1+DIAG_INDX-K, 1 )
	BNDEND(K)=MIN( (NUM_BNDS/ND)*(K-DIAG_INDX)+ND+DIAG_INDX-K,NUM_BNDS )
!
!	COUNT=COUNT+1			!Comment out for no diagnostics
	IF(COUNT .GT. 700)THEN
 	  PRINT_DIAGNOSTICS=.TRUE.
	ELSE
 	  PRINT_DIAGNOSTICS=.FALSE.
	END IF
	COMPUTE_LAM=.FALSE.
	FULL_ES=.TRUE.
	LUER=ERROR_LU()
	FREQ_INDX=0		!Not used
!
        IF(NNM .NE. 2 .AND. NNM .NE. 4)THEN
          WRITE(LUER,*)'Error NNM .NE. (2 .OR. 4) in SUB_CMF_LINE'
          STOP
        END IF
	IF(NP .NE. ND+NC-2)THEN
	  WRITE(LUER,*)'Error in SUB_CMF_LINE'
	  WRITE(LUER,*)'For this routine, NP must be ND+NC-2'
	  WRITE(LUER,*)'Advisable to use fine grid at outer boundary'
	  WRITE(LUER,*)'FG_COMP could be updated to handle ND+NC rays'
	  STOP
	END IF
C
C Increment emissivity due to electron scattering from the continuum.
C It is assumed that electron scattering from the line does not
C contribute significantly to the emissivity. The new ETA(I) is used
C by CMFJBAR
C
	DO I=1,ND
	  ETA(I)=ETA(I)+RJ(I)*CHI_SCAT(I)
	END DO
C
C NB - If COMPUTE_JEW is true, we set JEW to zero. It is then  initialized
C in the line subroutine (eg FORMSOL) to DNU*RJ.
C
	IF(COMPUTE_EW)THEN
	  IF(COMPUTE_JEW)THEN
	    CALL DP_ZERO(JEW,ND)
	  ELSE
	    READ(LU_JEW,REC=ACCESS_JEW)(JEW(I),I=1,ND),T1
	    IF(T1 .NE. FL)THEN
	      WRITE(LUER,*)'Error - incorrect reading of JEW'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP
	    END IF
	  END IF
	END IF
!
	IF(PRINT_DIAGNOSTICS)THEN
	  CALL WRITV(RJ,ND,'RJ',LUER)
	  CALL WRITV(ETA,ND,'ETA',LUER)
	  CALL WRITV(ETAL,ND,'ETAL',LUER)
	  CALL WRITV(CHI,ND,'CHI',LUER)
	  CALL WRITV(CHIL,ND,'CHIL',LUER)
	  CALL WRITV(CHI_SCAT,ND,'CHI_SCAT',LUER)
!
	  CALL WRITV(PF,NLF,'PF',LUER)
	  CALL WRITV(PROF,NLF,'PROF',LUER)
	  CALL WRITV(ERF,NLF,'ERF',LUER)
	  CALL WRITV(LFQW,NLF,'LFQW',LUER)
!
	  WRITE(LUER,*)'IC=',IC
	  WRITE(LUER,*)'FL=',FL
	  WRITE(LUER,*)'DBB=',DBB
	  WRITE(LUER,*)'DIF=',DIF
	  WRITE(LUER,*)'THK_CONT=',THK_CONT
	  WRITE(LUER,*)'THK_LINE=',THK_LINE
!
	END IF
!
C At present, higher accuracy computation for comoving frame
C calculation is not implemented. Will assume that can use
C accurate RJ from continuum source function. Thus eta and
C continuum source function are from accurate calculation.
C
	  IF(EDDINGTON)THEN
	    CALL TUNE(IONE,'MOMJBAR')
	    CALL DP_ZERO(JNU,ND*(NLF+1))
C
C As ETAL and CHIL are known, the FORMAL solution gives FEDD and GEDD
C directly. We only need to iterate to obtain a reliable EW value in
C the presence of electron scattering. Since the EW computation requires
C JEW to compute f, we need to save this for subsequent iterations.
C
	    DO I=1,ND
	      JNU(I,NLF+1)=JEW(I)
	    END DO
C
D	    WRITE(6,*)'Entering FG_COMP'
	    CALL TUNE(IONE,'FG_COMP')
	    CALL FG_COMP(ETA,CHI,CHI_SCAT,RJ,
	1                  CHIL,ETAL,V,SIGMA,R,P,
	1                  JNU,HNU,F_LINE,G_LINE,
	1                  AQW,HMIDQW,KQW,NMIDQW,
	1                  IN_HBC_LINE,HBC_LINE,NBC_LINE,
	1                  PF,PROF,LFQW,ERF,FL,
	1                  EW,CONT_INT,COMPUTE_EW,
	1                  DIF,DBB,IC,METHOD,
	1                  THK_CONT,THK_LINE,NLF,NC,NP,ND)
	    CALL TUNE(ITWO,'FG_COMP')
C
C We use TA for RADEQ, and TB for the FLUX vectors returned by the
C EW computation.
C
D	    WRITE(6,*)'Entering MOM_J_BAR'
	    CALL MOMJBAR(ETA,CHI,CHI_SCAT,THETA,RJ,CHIL,ETAL,
	1                  V,SIGMA,R,JBAR,ZNET,
	1                  JNU,HNU,F_LINE,G_LINE,
	1                  HBC_LINE,IN_HBC_LINE,NBC_LINE,JBLANK,HBLANK,
	1                  PF,PROF,LFQW,FL,DIF,DBB,IC,METHOD,
	1                  EW,CONT_INT,COMPUTE_EW,FULL_ES,
	1                  NLF,NC,NP,ND)
C
C Store J for EW computation in JEW for output to file.
C
	    DO I=1,ND
	      JEW(I)=JNU(I,NLF+1)
	    END DO
	    CALL TUNE(ITWO,'MOMJBAR')
	  ELSE
	    CALL TUNE(IONE,'FORMSOL')
	    CALL FORMSOL(ETA,CHI,CHI_SCAT,CHIL,ETAL,V,SIGMA,R,P,
	1                JBAR,ZNET,
	1                TA,TB,TC,COMPUTE_LAM,
	1                RJ,JBLANK,HBLANK,JEW,
	1                EW,CONT_INT,COMPUTE_EW,FULL_ES,
	1                AQW,HMIDQW,
	1                PF,PROF,LFQW,ERF,FL,DIF,DBB,IC,AMASS,
	1                THK_LINE,THK_CONT,NLF,NC,NP,ND,METHOD)
	    CALL TUNE(ITWO,'FORMSOL')
	  END IF
!
	IF(PRINT_DIAGNOSTICS)THEN
  	  CALL WRITV(JBAR,ND,'JBAR',LUER)
	END IF
C
C Store ZNET in ZNET_SIM for output on last iteration. We note simultaneous
C lines is not treated in CMF section.
C
	  DO K=1,ND
	    ZNET_SIM(K,1)=ZNET(K)
	  END DO
C
C
C Output JEW to file for use on subsequent iterations.
C
	  IF(COMPUTE_EW)THEN
	    WRITE(LU_JEW,REC=ACCESS_JEW)(JEW(I),I=1,ND),FL
	    ACCESS_JEW=ACCESS_JEW+1
	  END IF
C
C Increment the STEQ matrices due to Jmn line terms.
C
	  T1=FL*EMLIN
	  I=SIM_LINE_POINTER(1)
	  ID=VEC_ID(I)
	  MNL_F=VEC_MNL_F(I);     MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	  MNUP_F=VEC_MNUP_F(I);   MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
	  DO K=1,ND					!Equation depth
	    T2=ETAL_MAT(K,1)*ZNET(K)
	    SE(ID)%STEQ(MNUP,K)=SE(ID)%STEQ(MNUP,K) - T2/T1
	    SE(ID)%STEQ(MNL,K )=SE(ID)%STEQ(MNL,K) + T2/T1
	    STEQ_T(K)=STEQ_T(K) - T2
	  END DO
C
	  IF(COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION)THEN
D	    WRITE(6,*)'Entering LINEGEN'
	    CALL TUNE(IONE,'LINEGEN')
 	      INCLUDE 'LINEGEN.INC'
	    CALL TUNE(ITWO,'LINEGEN')
D	    WRITE(6,*)'Exiting LINEGEN'
	  ELSE IF(COMPUTE_BA)THEN
C
C If we are doing a lambda iteration, we are assuming that JBAR is fixed.
C Thus, dZ/dCHIL and dZ/dETAL is given by the following.
C
	    DO I=1,ND
	      VB(I)=-JBAR(I)/ETAL(I)
	      IF(NEG_OPACITY(I))VB(I)=0.0D0
	      VC(I)=JBAR(I)*CHIL(I)/ETAL(I)/ETAL(I)
	    END DO
C
C
C NB: We now include the temperature variation of L_STAR_RATIO and U_STAR_RATIO
C     in the variation of the net rate.
C
C Definitions:
C     RATE_FAC : Factor which multiplies the net rate to form the total
C                       rate from the upper level.
C     K          Depth index
C
	    NL=SIM_NL(1)
	    NUP=SIM_NUP(1)
	    I=SIM_LINE_POINTER(1)
	    ID=VEC_ID(I)
	    MNL_F=VEC_MNL_F(I)
	    MNUP_F=VEC_MNUP_F(I)
	    MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	    MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
	    MNT=SE(ID)%N_IV
!
	    T4=FL_SIM(1)*EMLIN
	    DO K=1,ND
	      L=GET_DIAG(K)
	      RATE_FAC=EINA(1)*POPS(NUP,K)*U_STAR_RATIO(K,1)
	      OPAC_FAC=OSCIL(1)*OPLIN
	      STIM_FAC=-GLDGU(1)*OPAC_FAC
	      EMIS_FAC=EINA(1)*FL_SIM(1)*EMLIN
	      dRATE_dUP=RATE_FAC*( ZNET(K)/POPS(NUP,K) +
	1               U_STAR_RATIO(K,1)*
	1              (STIM_FAC*VB(K)+EMIS_FAC*VC(K)) )
	      dRATE_dLOW=RATE_FAC*OPAC_FAC*VB(K)*
	1               L_STAR_RATIO(K,1)
	      dRATE_dT=RATE_FAC*
	1            (  OPAC_FAC*POPS(NL,K)*VB(K)*dL_RAT_dT(K,1)+
	1                POPS(NUP,K)*dU_RAT_dT(K,1)*(
	1                  EMIS_FAC*VC(K)+
	1                      STIM_FAC*VB(K) ) +
	1               ZNET(K)*dU_RAT_dT(K,1)
	1                          /U_STAR_RATIO(K,1)
	1            )
!
	      SE(ID)%BA(MNUP,MNUP,L,K)=SE(ID)%BA(MNUP,MNUP,L,K) - dRATE_dUP
	      SE(ID)%BA(MNUP,MNL,L,K) =SE(ID)%BA(MNUP,MNL,L,K) - dRATE_dLOW
	      SE(ID)%BA(MNUP,MNT,L,K) =SE(ID)%BA(MNUP,MNT,L,K) - dRATE_dT
	      SE(ID)%BA(MNL,MNUP,L,K) =SE(ID)%BA(MNL,MNUP,L,K)  + dRATE_dUP
	      SE(ID)%BA(MNL,MNL,L,K)  =SE(ID)%BA(MNL,MNL,L,K)  + dRATE_dLOW
	      SE(ID)%BA(MNL,MNT,L,K)  =SE(ID)%BA(MNL,MNT,L,K)  + dRATE_dT
!
	      BA_T(NL,L,K) =BA_T(NL,L,K)  - T4*dRATE_dLOW
	      BA_T(NUP,L,K)=BA_T(NUP,L,K) - T4*dRATE_dUP
	      BA_T(NT,L,K) =BA_T(NT,L,K)  - T4*dRATE_dT
	    END DO
C
	  END IF
!
D	  WRITE(6,*)'Exiting SUB_CMF_LINE'
C
C Write out rates and exit from NUP loop. Note that BA and STEQ have
C already been updated.
C
C	  GO TO 49000
C
	RETURN
	END
