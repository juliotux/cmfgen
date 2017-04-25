!
! Subroutine to compute the mean intensity. This routine replaces
! COMP_J_CONT.INC and was developed to faciltate inclusion of additional
! options without making CMFGEN_SUB continually larger.
!
! Available options for computing J: These should be the same as for DO_VAR_CONT.
! NB: DO_VAR CONT depends on results from this routine, however this routine is
! independent of DO_VAR_CONT.
!
!       CONT_VEL=.TRUE.			!Velocity field taken into account.
!            ACURATE=.TRUE.		!Enhanced spatial grid
!            ACCURATE=.FALSE.		!Regular spatial grid
!
!If the following options hold, we use Eddington factors and the enhanced grid.
!
!       CONT_VEL=.FALSE.  & THIS_FREQ_EXT=.TRUE.
!
!If the following options hold, we use Eddington factors and the regular grid.
!
!	CONT_VEL=.FALSE., EDDINGTON=.TRUE. & THIS_FREQ_EXT=.FALSE.
!
!If the following options hold we use ray-by ray solution for the full computaion.
!Spherical model, and no velocity terms.
!
!	CONT_VEL=.FALSE., EDDINGTON=.FALSE. & THIS_FREQ_EXT=.FALSE.
!
	SUBROUTINE COMP_J_BLANK(SECTION,EDDINGTON,FL,FREQ_INDX,FIRST_FREQ,LST_ITERATION,MAXCH,
	1                              LUER,LU_ES,LU_JCOMP,LU_EDD,ACCESS_F,
	1                              ND,NC,NP,NCF,NDEXT,NCEXT,NPEXT)
!
	USE ANG_QW_MOD
	USE MOD_CMFGEN
	USE OPAC_MOD
	USE RADIATION_MOD
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
! Altered : 17-Oct-2016 : H_CHK_OPTION added to moment routines.
!                            Moment routines are: MOM_J_CMF_V11.F, MOM_J_DDT_V4.F, and MOM_JREL_V8.F.
! Altered : 17-Feb-2015 : r^2.J and r^2.H now ouput on last iteration when using USE_LAM_ES option.
! Altered : 14-Dec-2014 : RSQHNU etc now set to one when not computing J. This avoids issues with
!                             possible NaNs.
! Altered : 16-Dec-2013 : CMF_FORM_SOL_V2 (non EXT option) no longer called when ND > 199.
!                             CMF_FORM_SOL_V2 is not parallelized and slows down large clumped models.
! Altered : 16-Feb-2006 : CMF_FORM_SOL_V2 used for last iteration when MAXCH<100, and
!                            not LAMBDA iteration. Sometimes it might be useful to
!                            change so that CMF_FORM_SOL_V2 is also called when LMABDA 
!                            iteration used. FG_COUNt was not being initialized.
! Finalized: 17-Dec-2004
!
	REAL*8 C_KMS
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	REAL*8 FL
	REAL*8 MAXCH
!
	INTEGER ACCESS_F
	INTEGER LU_EDD
	INTEGER LUER
	INTEGER LU_ES
	INTEGER LU_JCOMP
	INTEGER ND,NC,NP
	INTEGER NDEXT,NCEXT,NPEXT
	INTEGER NCF
	INTEGER FREQ_INDX
	CHARACTER(LEN=*) SECTION
!
	LOGICAL LST_ITERATION
	LOGICAL FIRST_FREQ
	LOGICAL EDDINGTON
!
! Local variables and arrays.
!
! Constants for opacity etc.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
	REAL*8, SAVE :: FL_OLD
	REAL*8 BNUE
	REAL*8 S1
	REAL*8 T1,T2
!
! FG_COUNT is used to determine the average number of calls to FG_J_CMF_V10 per
! frequency.
!
	INTEGER, SAVE :: FG_COUNT
	INTEGER I
	INTEGER J_IT_COUNTER
	LOGICAL NEW_FREQ
!
! See COMP_J_CONT.INC to see earlier changes to this routine.
!
! 
!
! Determine outer boundary confition. For CONT_VEL (i.e. full blanketing)
! THK_CONT is always set to RDRHK_CONT, and cannot change during a model run.
!
	  IF(ATM(1)%XzV_PRES)THEN		!Hydrogen
	    T1=ATM(1)%EDGEXzV_F(1)
	  ELSE IF(ATM(4)%XzV_PRES)THEN		!Helium II
	    T1=ATM(4)%EDGEXzV_F(2)
	  ELSE
	    T1=ATM(3)%EDGEXzV_F(5)		!Helium I
	  END IF
	  IF(RDTHK_CONT .AND. FL .GT. T1)THEN
	    THK_CONT=.TRUE.
	  ELSE
	    THK_CONT=.FALSE.
	  END IF
	  C_KMS=SPEED_OF_LIGHT(I)
!
! Option only avaialble for non-reativistic soultion, and no additinal
! points inserted.
!
	  OUT_BC_TYPE=1
	  IF(FL .LE. OUT_BC_PARAM_ONE)OUT_BC_TYPE=RD_OUT_BC_TYPE
!
       	IF(SECTION .EQ. 'CONTINUUM')THEN
	  CONT_VEL=.TRUE.
	  THK_CONT=RDTHK_CONT
	  IF(GLOBAL_LINE_SWITCH(1:5) .NE. 'BLANK' .AND. NO_VEL_FOR_CONTINUUM)THEN
	    CONT_VEL=.FALSE.
	  END IF
	ELSE
	  CONT_VEL=.FALSE.
	END IF
!
! Compute DBB and DDBBDT for diffusion approximation. DBB=dB/dR
! and DDBBDT= dB/dTR .
!
	T1=HDKT*FL/T(ND)
	T2=1.0D0-EMHNUKT(ND)
	BNUE=TWOHCSQ*( FL**3 )*EMHNUKT(ND)/T2
	DBB=TWOHCSQ*( FL**3 )*T1*DTDR/T(ND)*EMHNUKT(ND)/(T2**2)
	DDBBDT=DBB*(T1*(1.0D0+EMHNUKT(ND))/T2-2.0D0)/T(ND)
	HFLUX_AT_IB=DBB/CHI(ND)/3.0D0
!
! Switch to using CHI_CLUMP, ETA_CLUMP, and ESEC_CLUMP in case the model 
! has clumping.
!
	CHI_CLUMP(1:ND)=CHI(1:ND)*CLUMP_FAC(1:ND)
	ETA_CLUMP(1:ND)=ETA(1:ND)*CLUMP_FAC(1:ND)
	ESEC_CLUMP(1:ND)=ESEC(1:ND)*CLUMP_FAC(1:ND)
	CHI_SCAT_CLUMP(1:ND)=CHI_SCAT(1:ND)*CLUMP_FAC(1:ND)
!
! 
!                            
	IF(CONT_VEL .AND. USE_FIXED_J)THEN
	  CALL RD_CONT_J(FL,FREQ_INDX,FIRST_FREQ,LST_ITERATION,
	1          ACCURATE,LUER,LU_EDD,ACCESS_F,ND,NP)
	  RSQHNU=1.0D0; HFLUX_AT_OB=0.001D0; HFLUX_AT_IB=0.0D0
!
	ELSE IF(.NOT. CONT_VEL .AND. THIS_FREQ_EXT)THEN
C
C Solve for the mean intensity J . We can either solve for J with or without
C Eddington factors. Generally use Eddington factors when there is many 
C grid points.
C
	  CALL TUNE(IONE,'JFEAUEXT')
	  CALL EXTEND3OPAC(CHIEXT,ETAEXT,ESECEXT,COEF,INDX,
	1                      NDEXT,CHI_CLUMP,ETA_CLUMP,CHI_SCAT_CLUMP,ND)
C
	  DO I=1,NDEXT
	    ZETAEXT(I)=ETAEXT(I)/CHIEXT(I)
	    THETAEXT(I)=ESECEXT(I)/CHIEXT(I)
	  END DO
C
	  IF(SECTION .EQ. 'CONTINUUM' .AND. FREQ_INDX .EQ. 1)FEDD=0.0D0
	  IF(COMPUTE_EDDFAC)THEN
	    DO I=1,NDEXT
	      RJEXT(I)=0.0D0
	      FOLD(I)=0.0D0
	    END DO
	  ELSE
	    READ(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),T1
	    IF(ABS(T1/FL-1.0D0) .GT. 1.0D-10)THEN
	      WRITE(LUER,*)'Error - incorrect reading of EDDFACTOR in COMP_J_BLANK'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      WRITE(LUER,*)'You may need to delete EDDFACTOR'
	      STOP                
	    END IF
	  END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  J_IT_COUNTER=0
	  DO WHILE(INACCURATE)
	    DO I=1,NDEXT
	      SOURCEEXT(I)=ZETAEXT(I)+THETAEXT(I)*RJEXT(I)
	    END DO
	    S1=SOURCEEXT(1)
	    CALL FQCOMP_IBC(TA,TB,TC,XM,DTAU,REXT,Z,PEXT,QEXT,FEXT,
	1            SOURCEEXT,CHIEXT,dCHIdR,AQWEXT,KQWEXT,
	1            DBB,HBC_J,HBC_S,INBC,IC,
	1            THK_CONT,DIF,NCEXT,NDEXT,NPEXT,METHOD)
	    CALL JFEAU_IBC(TA,TB,TC,DTAU,REXT,RJEXT,QEXT,FEXT,
	1          ZETAEXT,THETAEXT,CHIEXT,DBB,IC,HBC_J,HBC_S,
	1          INBC,THK_CONT,DIF,NDEXT,METHOD)
C
C Update "inaccurate" iteration counter
C
	      J_IT_COUNTER=J_IT_COUNTER+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(J_IT_COUNTER .LT. 3 .OR. COMPUTE_EDDFAC)THEN	!Changed 8-Feb-95
	        T1=0.0D0
	        DO I=1,NDEXT
	          T1=MAX(ABS(FOLD(I)-FEXT(I)),T1)
	          FOLD(I)=FEXT(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
C
	      IF(J_IT_COUNTER .GT. 15)THEN
	         WRITE(LUER,*)'Possible error converging f - T1 is',T1
	         WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION 
	      	 INACCURATE=.FALSE.
	      END IF
	
	    END DO
C
C Put accurate calculation of J on old grid.
C
	    CALL UNGRID(RJ,ND,RJEXT,NDEXT,POS_IN_NEW_GRID)
	    CALL UNGRID(K_MOM,ND,FEXT,NDEXT,POS_IN_NEW_GRID)
	    DO I=1,ND
	      K_MOM(I)=K_MOM(I)*RJ(I)
	    END DO
C
C Optput Mean intensity for subsequent iterations.
C
	    WRITE(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),FL
C
C Update record for next frequency
	    ACCESS_F=ACCESS_F+1
C
	  CALL TUNE(ITWO,'JFEAUEXT')
C
C
C 
C
	ELSE IF(CONT_VEL .AND. ACCURATE)THEN
	  CALL TUNE(IONE,'JCONT_ACC')
C
C Interpolate the opacity and emissivity using a LINEAR interpolation
C law. CHIEXT etc. will contain the opacities etc. on the transfer grid
C with the clumping corrections. CHI_CLUMP etc refer to the appropriate
C quantities on the population grid.
C
	  CALL EXTEND3OPAC(CHIEXT,ETAEXT,ESECEXT,COEF,INDX,NDEXT,
	1              CHI_CLUMP,ETA_CLUMP,CHI_SCAT_CLUMP,ND)
C
C NB: CHI_PREV is used to refer to the continuum opacity at the previous
C frequency. Is does not need to be multiplied by CLUMP_FAC, as it is 
C compared directly to CHI_CONT. Since it is used for describing the
C variation in chi from one frequency to the next, we also do not need to
C use the extended vectors.
C
C For HBC and NBC only the first vector element is used.
C
	  CALL TUNE(IONE,'CONT_VEL')
	  NEW_FREQ=.TRUE.
	  IF(FIRST_FREQ)THEN
	    CHI_PREV(1:ND)=CHI(1:ND)
	    ETA_PREV(1:ND)=ETA(1:ND)
C
	    FEDD_PREV(1:NDEXT)=0.0D0		!Not required.
	    GEDD_PREV(1:NDEXT)=0.0D0
	    JNU_PREV(1:NDEXT)=0.0D0
	    N_ON_J_PREV(1:NDEXT)=0.0D0
	    RSQHNU_PREV(1:NDEXT)=0.0D0
C
	    HBC_PREV(:)=0.0D0		!1:3
	    NBC_PREV(:)=0.0D0		!1:3
	    HBC_CMF(:)=0.0D0		!1:3
	    NBC_CMF(:)=0.0D0		!1:3
	    FG_COUNT=0
	  ELSE
	    dLOG_NU=dLOG(FL_OLD/FL)
	    FEDD_PREV(1:NDEXT)=FEDD(1:NDEXT)
	    GEDD_PREV(1:NDEXT)=GEDD(1:NDEXT)
	    N_ON_J_PREV(1:NDEXT)=N_ON_J(1:NDEXT)
	    JNU_PREV(1:NDEXT)=RJEXT(1:NDEXT)
	    RSQHNU_PREV(1:NDEXT)=RSQHNU(1:NDEXT)
C
	    HBC_PREV(:)=HBC_CMF(:)
	    NBC_PREV(:)=NBC_CMF(:)
	  END IF
C
	  IF(COMPUTE_EDDFAC)THEN
	    IF(FIRST_FREQ)THEN
	      RJEXT(1:NDEXT)=0.0D0
	      RJEXT_ES(1:NDEXT)=0.0D0
	      FOLD(1:NDEXT)=0.0D0
	    END IF
	  ELSE
	    READ(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),T1
	    IF(ABS(T1/FL-1.0D0) .GT. 1.0D-10)THEN
	      WRITE(LUER,*)'Error - incorrect reading of EDDFACTOR in COMP_J_BLANK'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      WRITE(LUER,*)'You may need to delete EDDFACTOR'
	      STOP
	    END IF
	  END IF
C
C If we are using incoherent electron scattering, RJEXT_ES must be available.
C
	  IF(.NOT. COHERENT_ES)THEN
	    READ(LU_ES,REC=ACCESS_F)(RJEXT_ES(I),I=1,NDEXT),T1
	    IF(ABS(T1/FL-1.0D0) .GT. 1.0D-10)THEN
	      WRITE(LUER,*)'Error - incorrect reading of ES_J_CONV in COMO_J_BLANK'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      WRITE(LUER,*)'You may need to delete ES_J_CONV'
	      STOP
	    END IF
	  END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  J_IT_COUNTER=0
	  DO WHILE(INACCURATE)
C
	     IF(COHERENT_ES)THEN
	       TA(1:NDEXT)=ETAEXT(1:NDEXT) +
	1           ESECEXT(1:NDEXT)*RJEXT(1:NDEXT)
	     ELSE
	       TA(1:NDEXT)=ETAEXT(1:NDEXT) +
	1           ESECEXT(1:NDEXT)*RJEXT_ES(1:NDEXT)
	     END IF
C
C NB Using TA for ETA, TC for JNU_VEC, and TB for HNU_VEC
C
	     CALL TUNE(IONE,'FG_J_CMF_ACC')
	     CALL FG_J_CMF_V13(TA,CHIEXT,ESECEXT,
	1            VEXT,SIGMAEXT,REXT,PEXT,TC,FEDD,
	1            AQWEXT,HQWEXT,KQWEXT,NQWEXT,HMIDQWEXT,NMIDQWEXT,
	1            INBC,HBC_CMF(1),IPLUS,FL,dLOG_NU,
	1            INNER_BND_METH,DBB,IC,
	1            VDOP_VEC_EXT,DELV_FRAC_FG,REXT_FAC,
	1            METHOD,FG_SOL_OPTIONS,THK_CONT,
	1            FIRST_FREQ,NEW_FREQ,NCEXT,NPEXT,NDEXT)
	     CALL TUNE(ITWO,'FG_J_CMF_ACC')
	     FG_COUNT=FG_COUNT+1
C
	     IF(COHERENT_ES)THEN
	       TA(1:NDEXT)=ETAEXT(1:NDEXT)
	     ELSE
	       TA(1:NDEXT)=ETAEXT(1:NDEXT) +
	1            ESECEXT(1:NDEXT)*RJEXT_ES(1:NDEXT)
	     END IF
	     CALL TUNE(IONE,'MOM_J_CMF_ACC')
	     CALL MOM_J_CMF_V11(TA,CHIEXT,ESECEXT,VEXT,SIGMAEXT,REXT,
	1  	       RJEXT,RSQHNU,VDOP_VEC_EXT,DELV_FRAC_MOM,
	1              FL,dLOG_NU,INNER_BND_METH,DBB,IC,IB_STAB_FACTOR,
	1              N_TYPE,H_CHK_OPTION,METHOD,COHERENT_ES,OUT_BC_TYPE,
	1              FIRST_FREQ,NEW_FREQ,NCEXT,NPEXT,NDEXT)
	     CALL TUNE(ITWO,'MOM_J_CMF_ACC')
             IF(.NOT. DIF)HFLUX_AT_IB=0.5D0*IC*(0.5D0+INBC)-INBC*RJEXT(NDEXT)
             HFLUX_AT_OB=HBC_CMF(1)*RJEXT(1)
C
C We set NEW_FREQ to false so that FG_J_CMF continues to use the same
C AV_PREV and CV_PREV. NEW_FREQ must be set to true again outside the
C F iteration loop.
C
	     NEW_FREQ=.FALSE.
C
C Update "inaccurate" iteration counter
C
	      J_IT_COUNTER=J_IT_COUNTER+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(J_IT_COUNTER .LT. 20 .OR. COMPUTE_EDDFAC)THEN
	        T1=0.0D0
	        DO I=1,NDEXT
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)        
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
C
	      IF(J_IT_COUNTER .GT. 10)THEN
	         WRITE(LUER,*)'Possible error converging f - T1 is',T1
	         WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION 
	      	 INACCURATE=.FALSE.
	      END IF
	    END DO
C
C Output RJ for subsequent iterations.
C
	    WRITE(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),FL
C
C Store J on the normal mesh. No interpolation is involved here as
C it is assumed that the fine grid was created by the addition of extra
C points only.
C
	    CALL UNGRID(RJ,ND,RJEXT,NDEXT,POS_IN_NEW_GRID)
	    CALL UNGRID(RJ_ES,ND,RJEXT_ES,NDEXT,POS_IN_NEW_GRID)
	    CALL UNGRID(K_MOM,ND,FEDD,NDEXT,POS_IN_NEW_GRID)
C
C Compute K for use in computing mechanical energy loss.
C
	    K_MOM(1:ND)=RJ(1:ND)*K_MOM(1:ND)
C
	    IF(COHERENT_ES)THEN
	      SOURCE(1:ND)=ZETA(1:ND)+THETA(1:ND)*RJ(1:ND)
	    ELSE
	      SOURCE(1:ND)=ZETA(1:ND)+THETA(1:ND)*RJ_ES(1:ND)
	    END IF
C
C Update record for next frequency
C
	    ACCESS_F=ACCESS_F+1
	    FL_OLD=FL
!
! Note that TC is one the EXTENDED grid, hence we access its value at the
! inner boundary using NDEXT.
!
	  IF(LST_ITERATION)THEN
	    T1=ABS(RJ(1))+ABS(TC(1))
	    IF(T1 .NE. 0)T1=200.0D0*(RJ(1)-TC(1))/T1
	    T2=ABS(RJ(ND))+ABS(TC(NDEXT))
	    IF(T2 .NE. 0)T2=200.0D0*(RJ(ND)-TC(NDEXT))/T2
	    IF(FIRST_FREQ)THEN                 
	      OPEN(UNIT=LU_JCOMP,STATUS='UNKNOWN',FILE='J_COMP')
	      WRITE(LU_JCOMP,'(A)')' '
	      WRITE(LU_JCOMP,'(A)')'Comparison of J at Outer and Inner',
	1       ' boundaries computed using Moments and Ray techniques.'
	      WRITE(LU_JCOMP,'(A)')' '
	      WRITE(LU_JCOMP,
	1       '(3X,A,7X,A,7X,A,6X,A,5X,A,6X,A,5X,A,6X,A,5X,A)')
	1       'Indx','Nu','J(mom)','J(ray)','%Diff','HBC_CMF',
	1       'J(mom)','J(ray)','%Diff'
	    END IF
	    WRITE(LU_JCOMP,'(I7,ES16.6,2ES12.4,F10.2,3ES12.4,F10.2)')
	1                       FREQ_INDX,FL,
	1                       RJ(1),TC(1),T1,HBC_CMF(1),
	1                       RJ(ND),TC(NDEXT),T2
	    IF(FREQ_INDX .EQ. NCF)CLOSE(UNIT=LU_JCOMP)
	  END IF
	  IF(FREQ_INDX .EQ. NCF)THEN
	    WRITE(LUER,*)'Average number of calls to FG_J_CMF is',FLOAT(FG_COUNT)/FLOAT(NCF)
	    T1=ABS(RJ(1))+ABS(TC(1))
	    IF(T1 .NE. 0)T1=ABS(200.0D0*(RJ(1)-TC(1))/T1)
	    IF(T1 .GT. 100.0D0)THEN
	      WRITE(LUER,'(A)')'***************************************************************************'
	      WRITE(LUER,'(A)')' Error --- J(mom) and J(ray) differ by more than 100% for last frequency'
	      WRITE(LUER,'(A)')' It is STRONGLY suggested that you use a finer grid at the outer boundary'
	      WRITE(LUER,'(A)')' Tail J__COMP to see bundary error and/or use plt_jh'
	      WRITE(LUER,'(A)')'***************************************************************************'
	    ELSE IF(T1 .GT. 50.0D0)THEN
	      WRITE(LUER,'(A)')' Error --- J(mom) and J(ray) differ by more than 50% for last frequency'
	      WRITE(LUER,'(A)')' It is strongly suggested that you use a finer grid at the outer boundary'
	      WRITE(LUER,'(A)')' Tail J__COMP to see bundary error'
	    ELSE IF(T1 .GT. 20.0D0)THEN
	      WRITE(LUER,'(A)')' Error --- J(mom) and J(ray) differ by more than 20% for last frequency'
	      WRITE(LUER,'(A)')' Although this is nlikely to effect the colution, it is suggested that'
	      WRITE(LUER,'(A)')' you use a finer grid at the outer boundary'
	      WRITE(LUER,'(A)')' Tail J__COMP to see bundary error'
	    END IF
	  END IF
	  CALL TUNE(ITWO,'CONT_VEL')
C
C Set up for the compuation of the observes flux. LST_ITERATION is
C TRUE if FLUX_CAL_ONLY is true (single iteration with coherent,
C last iteration if non-coherent).
C
	    IF( (LST_ITERATION .AND. .NOT. LAMBDA_ITERATION .AND.
	1         MAXCH .LT. 100.0D0 .AND. .NOT. SN_MODEL) )THEN
C
C Quick and dirty method to ge an extended DENSITY vector. Will use TB in
C the call to CMF_FORM_SOL.
C
	       CALL EXTEND3OPAC(TA,TB,TC,COEF,INDX,NDEXT,
	1              DENSITY,DENSITY,DENSITY,ND)
C
	      IF(COHERENT_ES)THEN
	        TA(1:NDEXT)=ETAEXT(1:NDEXT)+ESECEXT(1:NDEXT)*RJEXT(1:NDEXT)
	      ELSE
	        TA(1:NDEXT)=ETAEXT(1:NDEXT)+ESECEXT(1:NDEXT)*RJEXT_ES(1:NDEXT)
	      END IF
C                             
C NB Using TA for ETA, U for P_OBS (temporay measure), I for NP_OBS.
C
	      CALL TUNE(IONE,'CMF_FORM_SOL')
	      CALL CMF_FORM_SOL_V2(TA,CHIEXT,ESECEXT,
	1                 TB,VEXT,SIGMAEXT,REXT,PEXT,
	1                 P_OBS,IPLUS,NP_OBS,NP_OBS_MAX,
	1                 MU_AT_RMAX,HQW_AT_RMAX,RMAX_OBS,V_AT_RMAX,
	1                 FL,dLOG_NU,DIF,DBB,IC,METHOD,
	1                 EXTEND_FRM_SOL,INSERT_FREQ_FRM_SOL,CMF_FORM_OPTIONS,
	1                 FRAC_DOP,V_DOP,dV_CMF_PROF,dV_CMF_WING,
	1                 FIRST_FREQ,NCEXT,NPEXT,NDEXT)
	      CALL TUNE(ITWO,'CMF_FORM_SOL')
	    ELSE IF(FIRST_FREQ)THEN
C
C So as defined for normal OBSFLUX calculation.
C
	      NP_OBS=NPEXT
	      P_OBS(1:NPEXT)=PEXT(1:NPEXT)
	      RMAX_OBS=R(1)
	      V_AT_RMAX=V(1)
	    END IF
!
	  CALL TUNE(ITWO,'JCONT_ACC')
C
C
C 
C
	ELSE IF(CONT_VEL)THEN
C
C NB: CHI_PREV is used to refer to the continuum opacity at the previous
C frequency. Is does not need to be multiplied by CLUMP_FAC, as it is compared
C directly to CHI_CONT.
C
C For HBC and NBC only the first vector element is used.
C
	  CALL TUNE(IONE,'CONT_VEL')
	  NEW_FREQ=.TRUE.
	  IF(FIRST_FREQ)THEN
	    DO I=1,ND
	      CHI_PREV(I)=CHI(I)
	      ETA_PREV(I)=ETA(I)
	      FEDD_PREV(I)=0.0D0		!Not required.
	      GEDD_PREV(I)=0.0D0
	      JNU_PREV(I)=0.0D0
	      N_ON_J_PREV(I)=0.0D0
	      RSQHNU_PREV(I)=0.0D0
	    END DO
	    HBC_PREV(:)=0.0D0		!1:3
	    NBC_PREV(:)=0.0D0		!1:3
	    HBC_CMF(:)=0.0D0		!1:3
	    NBC_CMF(:)=0.0D0		!1:3
	    FG_COUNT=0
	  ELSE
	    dLOG_NU=dLOG(FL_OLD/FL)
	    DO I=1,ND
	      FEDD_PREV(I)=FEDD(I)
	      GEDD_PREV(I)=GEDD(I)
	      N_ON_J_PREV(I)=N_ON_J(I)
	      JNU_PREV(I)=RJ(I)
	      RSQHNU_PREV(I)=RSQHNU(I)
	    END DO
	    HBC_PREV(:)=HBC_CMF(:)
	    NBC_PREV(:)=NBC_CMF(:)
	  END IF
C
	  IF(COMPUTE_EDDFAC)THEN
	    IF(FIRST_FREQ)THEN
	      DO I=1,ND
	        RJ(I)=0.0D0
	        RJ_ES(I)=0.0D0
	        FOLD(I)=0.0D0
	      END DO
	    END IF
	  ELSE
	    READ(LU_EDD,REC=ACCESS_F)(RJ(I),I=1,ND),T1
	    IF(ABS(T1/FL-1.0D0) .GT. 1.0D-10)THEN
	      WRITE(LUER,*)'Error - incorrect reading of EDDFACTOR in COMP_J_BLANK'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      WRITE(LUER,*)'You may need to delete EDDFACTOR'
	      STOP
	    END IF
	  END IF
C
C If we are using incoherent electron scattering, RJEXT_ES must be available.
C
	  IF(.NOT. COHERENT_ES)THEN
	    READ(LU_ES,REC=ACCESS_F)(RJ_ES(I),I=1,ND),T1
	    IF(T1 .NE. FL)THEN
	      WRITE(LUER,*)'Error - incorrect reading of ES_J_CONV in COMP_J_BLANK'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      WRITE(LUER,*)'You may need to delete ES_J_CONV'
	      STOP
	    END IF
	  END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  J_IT_COUNTER=0
	  DO WHILE(INACCURATE)
C
	     IF(COHERENT_ES)THEN
	       TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ(1:ND)
	     ELSE
	       TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ_ES(1:ND)
	     END IF
C
C NB Using TA for ETA, TC for JNU_VEC, and TB for HNU_VEC
C
	     CALL TUNE(IONE,'FG_J_CMF')
	     IF(PLANE_PARALLEL_NO_V)THEN
	        IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling FCOMP_PP'
	        SOURCE(1:ND)=TA(1:ND)/CHI_CLUMP(1:ND)
	        CALL FCOMP_PP_V2(R,TC,FEDD,SOURCE,CHI_CLUMP,IPLUS,HBC_CMF,
	1               NBC_CMF,INBC,DBB,IC,THK_CONT,DIF,ND,NC,METHOD)
	     ELSE IF(PLANE_PARALLEL)THEN
	        IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling PP_FORM_CMF_V2'
	        CALL PP_FORM_CMF_V2(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,
	1               TC,TB,FEDD,GEDD,N_ON_J,INBC,
	1               HBC_CMF(1),HBC_CMF(2),NBC_CMF(1),NBC_CMF(2),
	1               IPLUS,FL,dLOG_NU,DIF,DBB,IC,VDOP_VEC,DELV_FRAC_FG,
	1               METHOD,FG_SOL_OPTIONS,THK_CONT,INCL_INCID_RAD,
	1               FIRST_FREQ,NEW_FREQ,N_TYPE,NC,ND)
!
	     ELSE IF(USE_FORMAL_REL)THEN
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)
	1            WRITE(LUER,*)'Calling CMF_FORMAL_REL_V4 in COMP_J_BLANK'
	       CALL CMF_FORMAL_REL_V4
	1                 (TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,P,
	1                  TC,FEDD,HFLUX_AT_IB,HFLUX_AT_OB,IPLUS,
	1                  FL,dLOG_NU,BNUE,DBB,
	1                  INNER_BND_METH,THK_CONT,
	1                  VDOP_VEC,DELV_FRAC_FG,REXT_FAC,
	1                  METHOD,FIRST_FREQ,NEW_FREQ,NC,NP,ND)
!
	     ELSE 
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling FG_J_CMF_V13 in COMP_J_BLANK'
	       CALL FG_J_CMF_V13(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,P,
	1                  TC,FEDD,AQW,HQW,KQW,NQW,HMIDQW,NMIDQW,
	1                  INBC,HBC_CMF(1),IPLUS,FL,dLOG_NU,
	1                  INNER_BND_METH,DBB,IC,
	1                  VDOP_VEC,DELV_FRAC_FG,REXT_FAC,
	1                  METHOD,FG_SOL_OPTIONS,THK_CONT,
	1                  FIRST_FREQ,NEW_FREQ,NC,NP,ND)
!
	     END IF
	     CALL TUNE(ITWO,'FG_J_CMF')
	     FG_COUNT=FG_COUNT+1
C
	     IF(COHERENT_ES)THEN
	       TA(1:ND)=ETA_CLUMP(1:ND)
	     ELSE
	       TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ_ES(1:ND)
	     END IF
	     CALL TUNE(IONE,'MOM_J_CMF')
	     IF(PLANE_PARALLEL_NO_V)THEN
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling MOM_J_PP_V1'
	       CALL MOM_J_PP_V1(TA,CHI_CLUMP,CHI_SCAT_CLUMP,
	1                  R,FEDD,RJ,RSQHNU,HBC_CMF,NBC_CMF,INBC,
	1                  FL,DIF,DBB,IC,METHOD,COHERENT_ES,
	1                  IZERO,FIRST_FREQ,NEW_FREQ,ND)
	       HFLUX_AT_OB=HBC_CMF(1)*RJ(1)-HBC_CMF(2)
	       IF(.NOT. DIF)HFLUX_AT_OB=0.5D0*IC*(0.5D0+INBC)-INBC*RJ(ND)
	     ELSE IF(PLANE_PARALLEL)THEN
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling PP_MOM_CMF_V1'
	       CALL PP_MOM_CMF_V1(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,
	1                  FEDD,GEDD,N_ON_J,RJ,RSQHNU,
	1                  VDOP_VEC,DELV_FRAC_MOM,
	1                  INBC,HBC_CMF(1),HBC_CMF(2),NBC_CMF(1),NBC_CMF(2),
	1                  FL,dLOG_NU,DIF,DBB,IC,
	1                  N_TYPE,METHOD,COHERENT_ES,
	1                  FIRST_FREQ,NEW_FREQ,ND)
	       HFLUX_AT_OB=HBC_CMF(1)*RJ(1)-HBC_CMF(2)
	       IF(.NOT. DIF)HFLUX_AT_OB=0.5D0*IC*(0.5D0+INBC)-INBC*RJ(ND)
	     ELSE IF(USE_DJDT_RTE)THEN
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling MOM_J_DDT_V4'
	       CALL MOM_J_DDT_V4(TA,CHI_CLUMP,CHI_SCAT_CLUMP,
	1              V,R,FEDD,RJ,RSQHNU,DJDt_TERM,
	1              HFLUX_AT_IB,HFLUX_AT_OB,
	1              VDOP_VEC,DELV_FRAC_MOM,FL,dLOG_NU,DBB,
	1              H_CHK_OPTION,INNER_BND_METH,OUTER_BND_METH,
	1              METHOD,COHERENT_ES,FIRST_FREQ,NEW_FREQ,
	1              INCL_DJDT_TERMS,USE_DR4JDT,DJDT_RELAX_PARAM,NC,NP,ND,NCF)
	     ELSE IF(USE_J_REL)THEN
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling MOM_JREL_V8'
	       CALL MOM_JREL_V8(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,
	1             RJ,RSQHNU,HFLUX_AT_IB,HFLUX_AT_OB,
	1             VDOP_VEC,DELV_FRAC_MOM,
	1             FL,dLOG_NU,DBB,H_CHK_OPTION,IB_STAB_FACTOR,
	1             INNER_BND_METH,OUTER_BND_METH,
	1             METHOD,COHERENT_ES,N_TYPE,
	1             INCL_ADVEC_TERMS_IN_TRANS_EQ,INCL_REL_TERMS,FIRST_FREQ,ND)
	       IF(LST_ITERATION)THEN
	         DO I=1,ND
	           TA(I)=RJ(I)*R(I)*R(I)
	         END DO 
	         T1=HFLUX_AT_IB*R(ND)*R(ND)
	         T2=HFLUX_AT_OB/RJ(1)
	         CALL OUT_JH(TA,RSQHNU,T1,T2,FL,NCF,R,V,ND,FIRST_FREQ,'NORMAL')
	       END IF
	     ELSE IF(USE_LAM_ES)THEN
	       RJ(1:ND)=TC(1:ND)
	       IF(.NOT. USE_FORMAL_REL)THEN
	         IF(.NOT. DIF)HFLUX_AT_IB=0.5D0*IC*(0.5D0+INBC)-INBC*RJ(ND)
                 HFLUX_AT_OB=HBC_CMF(1)*RJ(1)
	       END IF
	       CALL GET_RSQH_REL(RSQHNU,R,V,FL,ND)
	       IF(LST_ITERATION)THEN
	         DO I=1,ND
	           TA(I)=RJ(I)*R(I)*R(I)
	         END DO 
	         T1=HFLUX_AT_IB*R(ND)*R(ND)
	         T2=HFLUX_AT_OB/RJ(1)
	         CALL OUT_JH(TA,RSQHNU,T1,T2,FL,NCF,R,V,ND,FIRST_FREQ,'NORMAL')
	       END IF
	     ELSE
	       IF(FIRST_FREQ .AND. J_IT_COUNTER .EQ. 0)WRITE(LUER,*)'Calling MOM_J_CMF_V11'
	       CALL MOM_J_CMF_V11(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,
	1              RJ,RSQHNU,VDOP_VEC,DELV_FRAC_MOM,
	1              FL,dLOG_NU,INNER_BND_METH,DBB,IC,IB_STAB_FACTOR,
	1              N_TYPE,H_CHK_OPTION,METHOD,COHERENT_ES,OUT_BC_TYPE,
	1              FIRST_FREQ,NEW_FREQ,NC,NP,ND)
	       IF(.NOT. DIF)HFLUX_AT_IB=0.5D0*IC*(0.5D0+INBC)-INBC*RJ(ND)
               HFLUX_AT_OB=HBC_CMF(1)*RJ(1)
	       IF(LST_ITERATION .AND. WRITE_JH)THEN
	         DO I=1,ND
	           TA(I)=RJ(I)*R(I)*R(I)
	         END DO
	         T1=HFLUX_AT_IB*R(ND)*R(ND)
	         T2=HFLUX_AT_OB/RJ(1)
	         CALL OUT_JH(TA,RSQHNU,T1,T2,FL,NCF,R,V,ND,FIRST_FREQ,'NORMAL')
	       END IF
	     END IF
	     CALL TUNE(ITWO,'MOM_J_CMF')
C
C We set NEW_FREQ to false so that FG_J_CMF continues to use the same
C AV_PREV and CV_PREV. NEW_FREQ must be set to true again outside the
C F iteration loop.
C
	     NEW_FREQ=.FALSE.
C
C Update "inaccurate" iteration counter
C
	      J_IT_COUNTER=J_IT_COUNTER+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(J_IT_COUNTER .LT. 20 .OR. COMPUTE_EDDFAC)THEN
	        T1=0.0D0
	        DO I=1,ND
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)        
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
C
	      IF(J_IT_COUNTER .GT. 20)THEN
	         WRITE(LUER,*)'Possible error converging f - T1 is',T1
	         WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION 
	      	 INACCURATE=.FALSE.
	      END IF
	    END DO
C
	    IF(COHERENT_ES)THEN
	      SOURCE(1:ND)=ZETA(1:ND)+THETA(1:ND)*RJ(1:ND)
	    ELSE
	      SOURCE(1:ND)=ZETA(1:ND)+THETA(1:ND)*RJ_ES(1:ND)
	    END IF
C
C Output RJ for subsequent iterations.
C
	    WRITE(LU_EDD,REC=ACCESS_F)(RJ(I),I=1,ND),FL
C
C Compute K for use in computing mechanical energy loss.
C
	    DO I=1,ND
	      K_MOM(I)=RJ(I)*FEDD(I)
	    END DO
C
C Update record for next frequency
C
	    ACCESS_F=ACCESS_F+1
	    FL_OLD=FL
C
C Set up for the compuation of the observes flux. LST_ITERATION is
C TRUE if FLUX_CAL_ONLY is true (single iteration with coherent,
C last iteration if non-coherent).
C
	    IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NO_V)THEN
!
! So as defined for normal OBSFLUX calculation. HQW_AT_RMAX is initially set 
! to JQW. Thus we need to multiply by MU to get the actual H weights at the
! outer boundary. For a plane-parallel atmosphere, RMAX_OBS is only scaling 
! constant. Setting its value to ND means that the observed luminosity should 
! correspond to the luminosity in VADAT (in absence of significant velocity 
! effects).
!
	      IF(FIRST_FREQ)THEN
	        NP_OBS=NC
	        CALL GAULEG(RZERO,RONE,MU_AT_RMAX,HQW_AT_RMAX,NC)
	        HQW_AT_RMAX(1:NC)=HQW_AT_RMAX(1:NC)*MU_AT_RMAX(1:NC)
	        RMAX_OBS=R(ND)
	        V_AT_RMAX=V(1)
	        IF(PLANE_PARALLEL_NO_V)V_AT_RMAX=0.0D0
	      END IF
	    ELSE IF( (LST_ITERATION .AND. .NOT. LAMBDA_ITERATION .AND.
	1         MAXCH .LT.  100.0D0 .AND. .NOT. SN_MODEL .AND. ND .LT. 200) )THEN
	      IF(COHERENT_ES)THEN
     	        TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ(1:ND)
	      ELSE
     	        TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ_ES(1:ND)
	      END IF
C                             
C NB Using TA for ETA, U for P_OBS (temporay measure), I for NP_OBS.
C
	      CALL TUNE(IONE,'CMF_FORM_SOL')
	      CALL CMF_FORM_SOL_V2(TA,CHI_CLUMP,CHI_SCAT_CLUMP,
	1                 DENSITY,V,SIGMA,R,P,
	1                 P_OBS,IPLUS,NP_OBS,NP_OBS_MAX,
	1                 MU_AT_RMAX,HQW_AT_RMAX,RMAX_OBS,V_AT_RMAX,
	1                 FL,dLOG_NU,DIF,DBB,IC,METHOD,
	1                 EXTEND_FRM_SOL,INSERT_FREQ_FRM_SOL,CMF_FORM_OPTIONS,
	1                 FRAC_DOP,V_DOP,dV_CMF_PROF,dV_CMF_WING,
	1                 FIRST_FREQ,NC,NP,ND)
	      CALL TUNE(ITWO,'CMF_FORM_SOL')
	    ELSE IF(FIRST_FREQ)THEN
!
! So as defined for normal OBSFLUX calculation. 
!
	      NP_OBS=NP
	      P_OBS(1:NP)=P(1:NP)
	      RMAX_OBS=R(1)
	      V_AT_RMAX=V(1)
	    END IF
C
	  IF(LST_ITERATION)THEN
	    T1=ABS(RJ(1))+ABS(TC(1))
	    IF(T1 .NE. 0)T1=200.0D0*(RJ(1)-TC(1))/T1
	    T2=ABS(RJ(ND))+ABS(TC(ND))
	    IF(T2 .NE. 0)T2=200.0D0*(RJ(ND)-TC(ND))/T2
	    IF(FIRST_FREQ)THEN                 
	      OPEN(UNIT=LU_JCOMP,STATUS='UNKNOWN',FILE='J_COMP')
	      WRITE(LU_JCOMP,'(A)')' '
	      WRITE(LU_JCOMP,'(A)')'Comparison of J at Outer and Inner',
	1       ' boundaries computed using Moments and Ray techniques.'
	      WRITE(LU_JCOMP,'(A)')' '
	      WRITE(LU_JCOMP,
	1       '(3X,A,14X,A,6X,A,6X,A,5X,A,5X,A,5X,A,7X,A,5X,A)')
	1       'Indx','Nu','J(mom)','J(ray)','%Diff','HBC_CMF',
	1       'J(mom)','J(ray)','%Diff'
	    END IF
	    WRITE(LU_JCOMP,'(I7,ES16.6,2ES12.4,F10.2,3ES12.4,F10.2)')
	1                       FREQ_INDX,FL,
	1                       RJ(1),TC(1),T1,HBC_CMF(1),
	1                       RJ(ND),TC(ND),T2
	    IF(FREQ_INDX .EQ. NCF)CLOSE(UNIT=LU_JCOMP)
	  END IF
	  IF(FREQ_INDX .EQ. NCF)THEN
	    WRITE(LUER,*)'Average number of calls to FG_J_CMF is',FLOAT(FG_COUNT)/FLOAT(NCF)
	    T1=ABS(RJ(1))+ABS(TC(1))
	    IF(T1 .NE. 0)T1=ABS(200.0D0*(RJ(1)-TC(1))/T1)
	    IF(T1 .GT. 100.0D0)THEN
	      WRITE(LUER,'(/,A)')' Error --- J(mom) and J(ray) differ by more than 1000% for last frequency'
	      WRITE(LUER,'(A)')' It is STRONGLY suggested that you use a finer grid at the outer boundary'
	      WRITE(LUER,'(A)')' Tail J__COMP to see bundary error and/or use plt_jh'
	    ELSE IF(T1 .GT. 50.0D0)THEN
	      WRITE(LUER,'(/,A)')' Error --- J(mom) and J(ray) differ by more than 50% for last frequency'
	      WRITE(LUER,'(A)')' It is strongly suggested that you use a finer grid at the outer boundary'
	      WRITE(LUER,'(A)')' Tail J__COMP to see bundary error'
	    ELSE IF(T1 .GT. 20.0D0)THEN
	      WRITE(LUER,'(/,A)')' Error --- J(mom) and J(ray) differ by more than 20% for last frequency'
	      WRITE(LUER,'(A)')' Although this is ulikely to effect the colution, it is suggested that'
	      WRITE(LUER,'(A)')' you use a finer grid at the outer boundary'
	      WRITE(LUER,'(A)')' Tail J__COMP to see bundary error'
	    END IF
	  END IF
	  CALL TUNE(ITWO,'CONT_VEL')
C
C 
C
	ELSE IF(EDDINGTON)THEN
C
C Calculation of "static" J in the continuum using Edington factors.
C
	  CALL TUNE(IONE,'JFEAU')
	  IF(SECTION .EQ. 'CONTINUUM' .AND. FREQ_INDX .EQ. 1)FEDD=0.0D0
	  IF(COMPUTE_EDDFAC)THEN
	    DO I=1,ND
	      RJ(I)=0.0D0
              FOLD(I)=FEDD(I)
	    END DO
	  ELSE
	    READ(LU_EDD,REC=ACCESS_F)(RJ(I),I=1,ND),T1
	    IF(T1 .NE. FL)THEN        
	      WRITE(LUER,'(/,A)')' Error - incorrect reading of EDDFACTOR in COMP_J_BLANK'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      WRITE(LUER,*)'You may need to delete EDDFACTOR'
	      STOP
	    END IF
	  END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  J_IT_COUNTER=0
	  DO WHILE(INACCURATE)
	    DO I=1,ND
	      SOURCE(I)=ZETA(I)+THETA(I)*RJ(I)
	    END DO
	    S1=SOURCE(1)
	    CALL FQCOMP_IBC(TA,TB,TC,XM,DTAU,R,Z,P,QEDD,FEDD,
	1            SOURCE,CHI_CLUMP,dCHIdR,AQW,KQW,DBB,HBC_J,HBC_S,
	1            INBC,IC,THK_CONT,DIF,NC,ND,NP,METHOD)
	    CALL JFEAU_IBC(TA,TB,TC,DTAU,R,RJ,QEDD,FEDD,
	1          ZETA,THETA,CHI_CLUMP,DBB,IC,HBC_J,HBC_S,
	1          INBC,THK_CONT,DIF,ND,METHOD)
C
C Update "inaccurate" iteration counter
C
	      J_IT_COUNTER=J_IT_COUNTER+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(J_IT_COUNTER .LT. 3 .OR. COMPUTE_EDDFAC)THEN
	        T1=0.0D0
	        DO I=1,ND
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF       
C
	      IF(J_IT_COUNTER .GT. 15)THEN
	         WRITE(LUER,*)'Possible error converging f - T1 is',T1
	         WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION 
	      	 INACCURATE=.FALSE.
	      END IF
	
	    END DO
C
	    DO I=1,ND
	      K_MOM(I)=RJ(I)*FEDD(I)
	    END DO
C
C Output mean intensity for subsequent iterations.
C
	    WRITE(LU_EDD,REC=ACCESS_F)(RJ(I),I=1,ND),FL
C
C Update record for next frequency
	    ACCESS_F=ACCESS_F+1
C
	  CALL TUNE(ITWO,'JFEAU')
	ELSE
C
C Calculation of "static" J in the continuum using Rybick method.
C
	  CALL TUNE(IONE,'JSOL')
	  CALL NEWJSOLD(TA,TB,TC,XM,WM,FB,RJ,DTAU,R,Z,P,
	1       ZETA,THETA,CHI_CLUMP,dCHIdR,AQW,
	1       THK_CONT,DIF,DBB,IC,NC,ND,NP,METHOD)
	  CALL TUNE(ITWO,'JSOL')
C
C Compute K_MOM.
C
	  DO I=1,ND
	    SOURCE(I)=ZETA(I)+THETA(I)*RJ(I)
	  END DO
	  S1=SOURCE(1)
	  CALL FQCOMP_IBC(TA,TB,TC,XM,DTAU,R,Z,P,QEDD,FEDD,
	1            SOURCE,CHI_CLUMP,dCHIdR,AQW,KQW,DBB,HBC_J,HBC_S,
	1            INBC,IC,THK_CONT,DIF,NC,ND,NP,METHOD)
	  DO I=1,ND
	    K_MOM(I)=K_MOM(I)*FEDD(I)
	  END DO
	END IF
!
	IF(ABS(RJ(1)) .GT. 1.0D+30)THEN
	  WRITE(LUER,*)' '
	  WRITE(LUER,*)' '
	  WRITE(LUER,*)'Error in comp_j_blank.f'
	  WRITE(LUER,*)'Mean intensity blowing up, which is due to an instabilty'
	  WRITE(LUER,*)'Try a finer grid at the outer boudary (preferred choice?)'
	  WRITE(LUER,*)'Alternatively try a different N_TYPE option'
	  STOP
	END IF
!
	RETURN
	END
