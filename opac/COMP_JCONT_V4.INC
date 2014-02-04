C
C Program section to compute the continuum intnesity. Opacities must have
C previuosly been evaluated.
C
C Altered 21-Apr-2001 - J [ (TC(NDEXT) ] was not being output to J_COMP correctly
C                         when extra points were being inserted. Now access
C                         TC(NDEXT), and use it before TC is corrupted.
C Altered 10-Apr-2000 - All methods for computing J now use CHI_CLUMP,
C                         ETA_CLUMP, and ESEC_CLUMP.
C Altered 11-Dec-1997 - Section for computing J in CONT_VEL section  with
C                         additional data points inserted. Designed to give
C                         higher accuracy, especially near ioization fronts.
C                         Additional points are inserted for ALL frequencies.
C                         For this option, ACCURATE must be set to TRUE.
C
C Altered 28-Jul-1997 - Correct DO_ LOOP variable now used to compute K_MOM
C                         in JSOL and JFEAU sections.
C Altered 17-May-1996 - Non-coherent electron scattering installed.
C Altered 16-Mar-1996 - FG_J_CMF and MOM_J_CMF updated to V5.
C
C 
C
C Determine outer boundary confition. For CONT_VEL (i.e. full blanketing)
C THK_CONT is always set to RDRHK_CONT, and cannot change during a model run.
C
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
C
       	IF(SECTION .EQ. 'CONTINUUM')THEN
	  CONT_VEL=.TRUE.
	  THK_CONT=RDTHK_CONT
	ELSE
	  CONT_VEL=.FALSE.
	END IF
C
C Compute DBB and DDBBDT for diffusion approximation. DBB=dB/dR
C and DDBBDT= dB/dTR .
C
	T1=HDKT*FL/T(ND)
	T2=1.0D0-EMHNUKT(ND)
	DBB=TWOHCSQ*( FL**3 )*T1*DTDR/T(ND)*EMHNUKT(ND)/(T2**2)
	DDBBDT=DBB*(T1*(1.0D0+EMHNUKT(ND))/T2-2.0D0)/T(ND)
C
C Switch to using CHI_CLUMP, ETA_CLUMP, and ESEC_CLUMP in case the model 
C has clumping.
C
	CHI_CLUMP(1:ND)=CHI(1:ND)*CLUMP_FAC(1:ND)
	ETA_CLUMP(1:ND)=ETA(1:ND)*CLUMP_FAC(1:ND)
	ESEC_CLUMP(1:ND)=ESEC(1:ND)*CLUMP_FAC(1:ND)
C
C 
C                            
	IF(.NOT. CONT_VEL .AND. THIS_FREQ_EXT)THEN
C
C Solve for the mean intensity J . We can either solve for J with or without
C Eddington factors. Generally use Eddington factors when there is many 
C grid points.
C
	  CALL TUNE(IONE,'JFEAUEXT')
	  CALL EXTEND3OPAC(CHIEXT,ETAEXT,ESECEXT,COEF,INDX,
	1                      NDEXT,CHI_CLUMP,ETA_CLUMP,ESEC_CLUMP,ND)
C
	  DO I=1,NDEXT
	    ZETAEXT(I)=ETAEXT(I)/CHIEXT(I)
	    THETAEXT(I)=ESECEXT(I)/CHIEXT(I)
	  END DO
C
	  IF(COMPUTE_EDDFAC)THEN
	    DO I=1,NDEXT
	      RJEXT(I)=0.0D0
	      FOLD(I)=0.0D0
	    END DO
	  ELSE
	    READ(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),T1
	    IF(T1 .NE. FL)THEN
	      WRITE(LUER,*)'Error - incorrect reading of mean intensity'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP                
	    END IF
	  END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  L=0
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
	      L=L+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(L .LT. 3 .OR. COMPUTE_EDDFAC)THEN	!Chnaged 8-Feb-95
	        T1=0.0D0
	        DO I=1,NDEXT
	          T1=MAX(ABS(FOLD(I)-FEXT(I)),T1)
	          FOLD(I)=FEXT(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
C
	      L=L+1
	      IF(L .GT. 15)THEN
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
	  CALL TUNE(IONE,'JCONT_ACC')
	  CALL EXTEND3OPAC(CHIEXT,ETAEXT,ESECEXT,COEF,INDX,NDEXT,
	1              CHI_CLUMP,ETA_CLUMP,ESEC_CLUMP,ND)
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
	    FG_COUNT=0.0D0
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
	    IF(T1 .NE. FL)THEN
	      WRITE(LUER,*)'Error - incorrect reading of'//
	1                  ' the mean intensity'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP
	    END IF
	  END IF
C
C If we are using incoherent electron scattering, RJEXT_ES must be available.
C
	  IF(.NOT. COHERENT_ES)THEN
	    READ(LU_ES,REC=ACCESS_F)(RJEXT_ES(I),I=1,NDEXT),T1
	    IF(T1 .NE. FL)THEN
	      WRITE(LUER,*)'Error - incorrect reading of'//
	1                ' the mean intensity'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP
	    END IF
	  END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  L=0
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
	     CALL TUNE(IONE,'FG_J_CMF')
	     CALL FG_J_CMF_V9(TA,CHIEXT,ESECEXT,
	1            VEXT,SIGMAEXT,REXT,PEXT,
	1            TC,TB,FEDD,GEDD,N_ON_J,
	1            AQWEXT,HMIDQWEXT,KQWEXT,NMIDQWEXT,
	1            INBC,HBC_CMF,NBC_CMF,IPLUS,
	1            FL,dLOG_NU,DIF,DBB,IC,
	1            METHOD,FG_SOL_OPTIONS,THK_CONT,
	1            FIRST_FREQ,NEW_FREQ,N_TYPE,
	1            NCEXT,NPEXT,NDEXT)
	     CALL TUNE(ITWO,'FG_J_CMF')
	     FG_COUNT=FG_COUNT+1
C
	     IF(COHERENT_ES)THEN
	       TA(1:NDEXT)=ETAEXT(1:NDEXT)
	     ELSE
	       TA(1:NDEXT)=ETAEXT(1:NDEXT) +
	1            ESECEXT(1:NDEXT)*RJEXT_ES(1:NDEXT)
	     END IF
	     CALL TUNE(IONE,'MOM_J_CMF')
	     CALL MOM_J_CMF_V5(TA,CHIEXT,ESECEXT,
	1              VEXT,SIGMAEXT,REXT,
	1  	       FEDD,GEDD,N_ON_J,
	1              FEDD_PREV,GEDD_PREV,N_ON_J_PREV,
	1              RJEXT,RSQHNU,JNU_PREV,RSQHNU_PREV,
	1              HBC_CMF,INBC,NBC_CMF,
	1              HBC_PREV,INBC_PREV,NBC_PREV,
	1              FL,dLOG_NU,DIF,DBB,IC,METHOD,COHERENT_ES,
	1              FIRST_FREQ,NEW_FREQ,
	1              NCEXT,NPEXT,NDEXT)
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
	      L=L+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(L .LT. 20 .OR. COMPUTE_EDDFAC)THEN
	        T1=0.0D0
	        DO I=1,NDEXT
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)        
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
C
	      IF(L .GT. 10)THEN
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
	    WRITE(LU_JCOMP,'(I7,1P3E12.4,0P,F10.2,1P,3E12.4,0P,F10.2)')
	1                       FREQ_INDX,FL,
	1                       RJ(1),TC(1),T1,HBC_CMF(1),
	1                       RJ(ND),TC(NDEXT),T2
	    IF(ML .EQ. NCF)CLOSE(UNIT=LU_JCOMP)
	  END IF
	  IF(ML .EQ. NCF)THEN
	     WRITE(LUER,*)'Average number of calls to FG_J_CMF is',FG_COUNT/NCF
	  END IF
C
C Set up for the compuation of the observes flux. LST_ITERATION is
C TRUE if FLUX_CAL_ONLY is true (single iteration with coherent,
C last iteration if non-coherent).
C
	    IF( (LST_ITERATION .AND. .NOT. LAMBDA_ITERATION .AND.
	1         MAXCH .LT. VAL_DO_NG) )THEN
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
	      CALL CMF_FORM_SOL(TA,CHIEXT,ESECEXT,
	1                 TB,VEXT,SIGMAEXT,REXT,PEXT,
	1                 P_OBS,IPLUS,NP_OBS,NP_OBS_MAX,
	1                 MU_AT_RMAX,HQW_AT_RMAX,RMAX_OBS,V_AT_RMAX,
	1                 FL,dLOG_NU,DIF,DBB,IC,METHOD,
	1                 EXTEND_FRM_SOL,INSERT_FREQ_FRM_SOL,
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
	  CALL TUNE(ITWO,'CONT_VEL')
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
	    FG_COUNT=0.0D0
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
	    IF(T1 .NE. FL)THEN
	      WRITE(LUER,*)'Error - incorrect reading of'//
	1                  ' the mean intensity'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP
	    END IF
	  END IF
C
C If we are using incoherent electron scattering, RJEXT_ES must be available.
C
	  IF(.NOT. COHERENT_ES)THEN
	    READ(LU_ES,REC=ACCESS_F)(RJ_ES(I),I=1,ND),T1
	    IF(T1 .NE. FL)THEN
	      WRITE(LUER,*)'Error - incorrect reading of'//
	1                ' the mean intensity'
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP
	    END IF
	  END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  L=0
	  DO WHILE(INACCURATE)
C
	     IF(COHERENT_ES)THEN
	       TA(1:ND)=ETA_CLUMP(1:ND)+ESEC_CLUMP(1:ND)*RJ(1:ND)
	     ELSE
	       TA(1:ND)=ETA_CLUMP(1:ND)+ESEC_CLUMP(1:ND)*RJ_ES(1:ND)
	     END IF
C
C NB Using TA for ETA, TC for JNU_VEC, and TB for HNU_VEC
C
	     CALL TUNE(IONE,'FG_J_CMF')
	     CALL FG_J_CMF_V9(TA,CHI_CLUMP,ESEC_CLUMP,V,SIGMA,R,P,
	1                  TC,TB,FEDD,GEDD,N_ON_J,
	1                  AQW,HMIDQW,KQW,NMIDQW,
	1                  INBC,HBC_CMF,NBC_CMF,
	1                  IPLUS,FL,dLOG_NU,DIF,DBB,IC,
	1                  METHOD,FG_SOL_OPTIONS,THK_CONT,
	1                  FIRST_FREQ,NEW_FREQ,N_TYPE,NC,NP,ND)
	     CALL TUNE(ITWO,'FG_J_CMF')
	     FG_COUNT=FG_COUNT+1
C
	     IF(COHERENT_ES)THEN
	       TA(1:ND)=ETA_CLUMP(1:ND)
	     ELSE
	       TA(1:ND)=ETA_CLUMP(1:ND)+ESEC_CLUMP(1:ND)*RJ_ES(1:ND)
	     END IF
	     CALL TUNE(IONE,'MOM_J_CMF')
	     CALL MOM_J_CMF_V5(TA,CHI_CLUMP,ESEC_CLUMP,V,SIGMA,R,
	1  	       FEDD,GEDD,N_ON_J,
	1              FEDD_PREV,GEDD_PREV,N_ON_J_PREV,
	1              RJ,RSQHNU,JNU_PREV,RSQHNU_PREV,
	1              HBC_CMF,INBC,NBC_CMF,
	1              HBC_PREV,INBC_PREV,NBC_PREV,
	1              FL,dLOG_NU,DIF,DBB,IC,METHOD,COHERENT_ES,
	1              FIRST_FREQ,NEW_FREQ,NC,NP,ND)
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
	      L=L+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(L .LT. 20 .OR. COMPUTE_EDDFAC)THEN
	        T1=0.0D0
	        DO I=1,ND
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)        
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
C
	      IF(L .GT. 10)THEN
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
	    IF( (LST_ITERATION .AND. .NOT. LAMBDA_ITERATION .AND.
	1         MAXCH .LT. VAL_DO_NG) )THEN
	      IF(COHERENT_ES)THEN
     	        TA(1:ND)=ETA_CLUMP(1:ND)+ESEC_CLUMP(1:ND)*RJ(1:ND)
	      ELSE
     	        TA(1:ND)=ETA_CLUMP(1:ND)+ESEC_CLUMP(1:ND)*RJ_ES(1:ND)
	      END IF
C                             
C NB Using TA for ETA, U for P_OBS (temporay measure), I for NP_OBS.
C
	      CALL TUNE(IONE,'CMF_FORM_SOL')
	      CALL CMF_FORM_SOL(TA,CHI_CLUMP,ESEC_CLUMP,
	1                 DENSITY,V,SIGMA,R,P,
	1                 P_OBS,IPLUS,NP_OBS,NP_OBS_MAX,
	1                 MU_AT_RMAX,HQW_AT_RMAX,RMAX_OBS,V_AT_RMAX,
	1                 FL,dLOG_NU,DIF,DBB,IC,METHOD,
	1                 EXTEND_FRM_SOL,INSERT_FREQ_FRM_SOL,
	1                 FRAC_DOP,V_DOP,dV_CMF_PROF,dV_CMF_WING,
	1                 FIRST_FREQ,NC,NP,ND)
	      CALL TUNE(ITWO,'CMF_FORM_SOL')
	    ELSE IF(FIRST_FREQ)THEN
C
C So as defined for normal OBSFLUX calculation.
C
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
	1       '(3X,A,7X,A,7X,A,6X,A,5X,A,6X,A,5X,A,6X,A,5X,A)')
	1       'Indx','Nu','J(mom)','J(ray)','%Diff','HBC_CMF',
	1       'J(mom)','J(ray)','%Diff'
	    END IF
	    WRITE(LU_JCOMP,'(I7,1P3E12.4,0P,F10.2,1P,3E12.4,0P,F10.2)')
	1                       FREQ_INDX,FL,
	1                       RJ(1),TC(1),T1,HBC_CMF(1),
	1                       RJ(ND),TC(ND),T2
	    IF(ML .EQ. NCF)CLOSE(UNIT=LU_JCOMP)
	  END IF
	  IF(ML .EQ. NCF)THEN
	     WRITE(LUER,*)'Average number of calls to FG_J_CMF is',FG_COUNT/NCF
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
	  IF(COMPUTE_EDDFAC)THEN
	    DO I=1,ND
	      RJ(I)=0.0D0
              FOLD(I)=FEDD(I)
	    END DO
	  ELSE
	    READ(LU_EDD,REC=ACCESS_F)(RJ(I),I=1,ND),T1
	    IF(T1 .NE. FL)THEN        
	      WRITE(LUER,*)'Error - incorrect reading of Mean Intensity' 
	      WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP
	    END IF
	  END IF
C
C We will do this twice, so that F is of higher accuracy.
C
	  INACCURATE=.TRUE.
	  L=0
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
	      L=L+1
C
C Check if F has converged.
C
	      INACCURATE=.FALSE.
	      IF(L .LT. 3 .OR. COMPUTE_EDDFAC)THEN
	        T1=0.0D0
	        DO I=1,ND
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF       
C
	      L=L+1
	      IF(L .GT. 15)THEN
	         WRITE(LUER,*)'Possible error converging f - T1 is',T1
	         WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION 
	      	 INACCURATE=.FALSE.
	      END IF
	
	    END DO
C
	    DO I=1,ND
	      K_MOM(I)=K_MOM(I)*FEDD(I)
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
C
