!
! Program section to compute the continuum intensity. Opacities must have
! previuosly been evaluated.
!
! Altered 22-Apr-2015 - Added USE_DR4JDT option. Changed to MOM_J_DDT_V3
! Altered 09-Jun-2015 - The frequency in EDDFACTOR etc now only have to be accurate to 1 part in 1.0D-14.
!                          (after Ralph Sutherland, although % comparison.)
! Altered 14-mar-2013 - When useing EDDINGTON we now call jfeau_ibc_v2.f and fqcomp_ibc_v2.f.
!                          These handle a hollow core/zero flux inner boundary condition.
! Altered 10-Apr-2011 - Bug fix: HFLUX_AT_OB and HFLUX_AT_IB were not
!                         being correctly set for normal wind calculation. Only
!                         effects  output to FLUX_FILE.
! Altered 24-Jan-2009 - Changes to call MOM_JREL_V6
!                       Call to CMF_FORMAL_REL_V3
!                       Variable HFLUX_AT_OB, HFLUX_AT_IB introduced.
! Altered 31-Jul-2007 - Call to MOM_JREL_V3 etc installed.
! Altered 13-Apr-2007 - Fixed routine so that V_AT_RMAX=0 for plane-parallel model
!                          without V, even if V is set. Useful for testing.
! Altered 20-Mar-2007 - Installed plane-parallel routines for Accurate grid.
! Altered 26-Apr-2006 - Installed plane-parallel routines (no EXT only)
! Altered 21-Apr-2001 - J [ (TC(NDEXT) ] was not being output to J_COMP correctly
!                         when extra points were being inserted. Now access
!                         TC(NDEXT), and use it before TC is corrupted.
! Altered 10-Apr-2000 - All methods for computing J now use CHI_CLUMP,
!                         ETA_CLUMP, and ESEC_CLUMP.
! Altered 11-Dec-1997 - Section for computing J in CONT_VEL section  with
!                         additional data points inserted. Designed to give
!                         higher accuracy, especially near ioization fronts.
!                         Additional points are inserted for ALL frequencies.
!                         For this option, ACCURATE must be set to TRUE.
!
! Altered 28-Jul-1997 - Correct DO_ LOOP variable now used to compute K_MOM
!                         in JSOL and JFEAU sections.
! Altered 17-May-1996 - Non-coherent electron scattering installed.
! Altered 16-Mar-1996 - FG_J_CMF and MOM_J_CMF updated to V5.
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
!
       	IF(SECTION .EQ. 'CONTINUUM')THEN
	  CONT_VEL=.TRUE.
	  THK_CONT=RDTHK_CONT
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
	IF(SN_MODEL .AND. .NOT. DIF)THEN
	  BNUE=0.0D0; DBB=0.0D0; DDBBDT=0.0D0
	  IF(FIRST_FREQ)WRITE(LUER,*)'Setting flux for hollow core in CMF_FLUX_SUB_V5'
	END IF
!
! Switch to using CHI_CLUMP, ETA_CLUMP, and ESEC_CLUMP in case the model
! has clumping.
!
	CHI_CLUMP(1:ND)=CHI(1:ND)*CLUMP_FAC(1:ND)
	ETA_CLUMP(1:ND)=ETA(1:ND)*CLUMP_FAC(1:ND)
	ESEC_CLUMP(1:ND)=ESEC(1:ND)*CLUMP_FAC(1:ND)
	CHI_RAY_CLUMP(1:ND)=CHI_RAY(1:ND)*CLUMP_FAC(1:ND)
	CHI_SCAT_CLUMP(1:ND)=CHI_SCAT(1:ND)*CLUMP_FAC(1:ND)
!
! 
!
	IF(CONT_VEL .AND. USE_FIXED_J)THEN
	  IF(ACCURATE)THEN
	    CALL RD_CONT_J_OBS(RJEXT,FL,FREQ_INDX,FIRST_FREQ,
	1          LST_ITERATION,LUER,LU_EDD,ACCESS_F,NDEXT)
	    CALL UNGRID(RJ,ND,RJEXT,NDEXT,POS_IN_NEW_GRID)
	    RSQHNU(1:NDEXT)=1.0D-20
	  ELSE
	    CALL RD_CONT_J_OBS(RJ,FL,FREQ_INDX,FIRST_FREQ,
	1          LST_ITERATION,LUER,LU_EDD,ACCESS_F,ND)
	    RSQHNU(1:ND)=1.0D-20
	  END IF
!
! Evaluate other quantities that are set but not really needed.
!
	  RJ_ES(1:ND)=RJ(1:ND)
	  K_MOM(1:ND)=RJ(1:ND)/3.0D0
          SOURCE(1:ND)=ZETA(1:ND)+THETA(1:ND)*RJ(1:ND)
!
! So as defined for normal OBSFLUX calculation.
!
	  IF(FIRST_FREQ)THEN
	    NP_OBS=NP
	    P_OBS(1:NP)=P(1:NP)
	    RMAX_OBS=R(1)
	    V_AT_RMAX=V(1)
	  END IF
	  IPLUS(1:NP)=1.0D-10             !Arbitrary value
	  HBC_CMF(1)=1.0D0
!
	ELSE IF(.NOT. CONT_VEL .AND. THIS_FREQ_EXT)THEN
!
! Solve for the mean intensity J . We can either solve for J with or without
! Eddington factors. Generally use Eddington factors when there is many
! grid points.
!
	  CALL TUNE(IONE,'JFEAUEXT')
	  CALL EXTEND3OPAC(CHIEXT,ETAEXT,ESECEXT,COEF,INDX,
	1                      NDEXT,CHI_CLUMP,ETA_CLUMP,ESEC_CLUMP,ND)
	  CALL EXTEND3OPAC(CHIEXT,CHI_RAY_EXT,CHI_SCAT_EXT,COEF,INDX,
	1                      NDEXT,CHI_CLUMP,CHI_RAY_CLUMP,CHI_SCAT_CLUMP,ND)
!
	  DO I=1,NDEXT
	    ZETAEXT(I)=ETAEXT(I)/CHIEXT(I)
	    THETAEXT(I)=CHI_SCAT_EXT(I)/CHIEXT(I)
	  END DO
!
	  IF(COMPUTE_EDDFAC)THEN
	    DO I=1,NDEXT
	      RJEXT(I)=0.0D0
	      FOLD(I)=0.0D0
	    END DO
	  ELSE
	    READ(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),T1
	    IF(ABS(T1-FL)/FL .GT. 1.0D-14)THEN
	      WRITE(LUER,*)'Error - 1 incorrect reading of mean intensity'
	      WRITE(LUER,*)'Intensity is ',FL,'Old Intensity is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP
	    END IF
	  END IF
!
! We will do this twice, so that F is of higher accuracy.
!
	  INACCURATE=.TRUE.
	  L=0
	  DO WHILE(INACCURATE)
	    DO I=1,NDEXT
	      SOURCEEXT(I)=ZETAEXT(I)+THETAEXT(I)*RJEXT(I)
	    END DO
	    S1=SOURCEEXT(1)
	    CALL FQCOMP_IBC_V2(TA,TB,TC,XM,DTAU,REXT,Z,PEXT,QEXT,FEXT,
	1            SOURCEEXT,CHIEXT,dCHIdR,AQWEXT,KQWEXT,
	1            DBB,HBC_J,HBC_S,INBC,IC,
	1            THK_CONT,INNER_BND_METH,NCEXT,NDEXT,NPEXT,METHOD)
	    CALL JFEAU_IBC_V2(TA,TB,TC,DTAU,REXT,RJEXT,QEXT,FEXT,
	1          ZETAEXT,THETAEXT,CHIEXT,DBB,IC,HBC_J,HBC_S,
	1          INBC,THK_CONT,INNER_BND_METH,NDEXT,METHOD)
!
! Update "inaccurate" iteration counter
!
	      L=L+1
!
! Check if F has converged.
!
	      INACCURATE=.FALSE.
	      IF(L .LT. 3 .OR. COMPUTE_EDDFAC)THEN	!Chnaged 8-Feb-95
	        T1=0.0D0
	        DO I=1,NDEXT
	          T1=MAX(ABS(FOLD(I)-FEXT(I)),T1)
	          FOLD(I)=FEXT(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
!
	      L=L+1
	      IF(L .GT. 15)THEN
	         WRITE(LUER,*)'Possible error converging f - T1 is',T1
	         WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION
	      	 INACCURATE=.FALSE.
	      END IF

	    END DO
!
! Put accurate calculation of J on old grid.
!
	    CALL UNGRID(RJ,ND,RJEXT,NDEXT,POS_IN_NEW_GRID)
	    CALL UNGRID(K_MOM,ND,FEXT,NDEXT,POS_IN_NEW_GRID)
	    DO I=1,ND
	      K_MOM(I)=K_MOM(I)*RJ(I)
	    END DO
!
! Optput Mean intensity for subsequent iterations.
!
	    WRITE(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),FL
!
! Update record for next frequency
	    ACCESS_F=ACCESS_F+1
!
	  CALL TUNE(ITWO,'JFEAUEXT')
!
! 
!
	ELSE IF(CONT_VEL .AND. ACCURATE)THEN
!
! Interpolate the opacity and emissivity using a LINEAR interpolation
! law. CHIEXT etc. will contain the opacities etc. on the transfer grid
! with the clumping corrections. CHI_CLUMP etc refer to the appropriate
! quantities on the population grid.
!
	  CALL EXTEND3OPAC(CHIEXT,ETAEXT,ESECEXT,COEF,INDX,NDEXT,
	1              CHI_CLUMP,ETA_CLUMP,ESEC_CLUMP,ND)
	  CALL EXTEND3OPAC(CHIEXT,CHI_RAY_EXT,CHI_SCAT_EXT,COEF,INDX,
	1                      NDEXT,CHI_CLUMP,CHI_RAY_CLUMP,CHI_SCAT_CLUMP,ND)
!
! NB: CHI_PREV is used to refer to the continuum opacity at the previous
! frequency. Is does not need to be multiplied by CLUMP_FAC, as it is
! compared directly to CHI_CONT. Since it is used for describing the
! variation in chi from one frequency to the next, we also do not need to
! use the extended vectors.
!
! For HBC and NBC only the first vector element is used.
!
	  CALL TUNE(IONE,'CONT_VEL')
	  NEW_FREQ=.TRUE.
	  IF(FIRST_FREQ)THEN
	    CHI_PREV(1:ND)=CHI(1:ND)
	    ETA_PREV(1:ND)=ETA(1:ND)
!
	    FEDD_PREV(1:NDEXT)=0.0D0		!Not required.
	    GEDD_PREV(1:NDEXT)=0.0D0
	    JNU_PREV(1:NDEXT)=0.0D0
	    N_ON_J_PREV(1:NDEXT)=0.0D0
	    RSQHNU_PREV(1:NDEXT)=0.0D0
!
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
!
	    HBC_PREV(:)=HBC_CMF(:)
	    NBC_PREV(:)=NBC_CMF(:)
	  END IF
!
	  IF(COMPUTE_EDDFAC)THEN
	    IF(FIRST_FREQ)THEN
	      RJEXT(1:NDEXT)=0.0D0
	      RJEXT_ES(1:NDEXT)=0.0D0
	      FOLD(1:NDEXT)=0.0D0
	    END IF
	  ELSE
	    READ(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),T1
	    IF(ABS(T1-FL)/FL .GT. 1.0D-14)THEN
	      WRITE(LUER,*)'Error - 2 incorrect reading of'//
	1                  ' the mean intensity'
	      WRITE(LUER,*)'Intensity is ',FL,'Old Intensity is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP
	    END IF
	  END IF
!
! If we are using incoherent electron scattering, RJEXT_ES must be available.
!
	  IF(.NOT. COHERENT_ES)THEN
	    READ(LU_ES,REC=ACCESS_F)(RJEXT_ES(I),I=1,NDEXT),T1
	    IF(ABS(T1-FL)/FL .GT. 1.0D-14)THEN
	      WRITE(LUER,*)'Error - 3 incorrect reading of'//
	1                ' the mean intensity'
	      WRITE(LUER,*)'Intensity is ',FL,'Old Intensity is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP
	    END IF
	  END IF
!
! We will do this twice, so that F is of higher accuracy.
!
	  INACCURATE=.TRUE.
	  L=0
	  DO WHILE(INACCURATE)
!
! We assume the Rayleigh scattering to be coherent in the CMF.
!
	     IF(COHERENT_ES)THEN
	       TA(1:NDEXT)=ETAEXT(1:NDEXT) +
	1           CHI_SCAT_EXT(1:NDEXT)*RJEXT(1:NDEXT)
	     ELSE
	       TA(1:NDEXT)=ETAEXT(1:NDEXT) +
	1           ESECEXT(1:NDEXT)*RJEXT_ES(1:NDEXT) +
	1           CHI_RAY_EXT(1:NDEXT)*RJEXT(1:NDEXT)
	     END IF
!
! NB Using TA for ETA, TC for JNU_VEC, and TB for HNU_VEC
!
	     FG_COUNT=FG_COUNT+1
	     CALL TUNE(IONE,'FG_J_CMF_ACC')
	     IF(PLANE_PARALLEL_NO_V)THEN
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling FCOMP_PP_V2 with accurate grid'
	       SOURCEEXT(1:NDEXT)=TA(1:NDEXT)/CHIEXT(1:NDEXT)
	       CALL FCOMP_PP_V2(REXT,TC,FEDD,SOURCEEXT,CHIEXT,IPLUS,HBC_CMF,
	1               NBC_CMF,INBC,DBB,IC,THK_CONT,DIF,NDEXT,NC,METHOD)
	     ELSE IF(PLANE_PARALLEL)THEN
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling PP_FORM_CMF_V2 with accurate grid'
	       CALL PP_FORM_CMF_V2(TA,CHIEXT,CHI_SCAT_EXT,VEXT,SIGMAEXT,REXT,
	1               TC,TB,FEDD,GEDD,N_ON_J,INBC,
	1               HBC_CMF(1),HBC_CMF(2),NBC_CMF(1),NBC_CMF(2),
	1               IPLUS,FL,dLOG_NU,DIF,DBB,IC,VDOP_VEC_EXT,DELV_FRAC_FG,
	1               METHOD,FG_SOL_OPTIONS,THK_CONT,L_FALSE,
	1               FIRST_FREQ,NEW_FREQ,N_TYPE,NC,NDEXT)
             ELSE IF(USE_FORMAL_REL)THEN
               IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling CMF_FORMAL_REL_V4 with accurate grid in COMP_JCONT_V4.INC'
               CALL CMF_FORMAL_REL_V4
	1              (TA,CHIEXT,CHI_SCAT_EXT,VEXT,SIGMAEXT,REXT,PEXT,
	1              TC,FEDD,HFLUX_AT_IB,HFLUX_AT_OB,IPLUS,FL,dLOG_NU,
	1              BNUE,DBB,INNER_BND_METH,THK_CONT,
	1              VDOP_VEC_EXT,DELV_FRAC_FG,REXT_FAC,
	1              METHOD,FIRST_FREQ,NEW_FREQ,NCEXT,NPEXT,NDEXT)
	     ELSE
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling FG_J_CMF_V12 with accurate grid'
	       CALL FG_J_CMF_V12(TA,CHIEXT,CHI_SCAT_EXT,
	1              VEXT,SIGMAEXT,REXT,PEXT,TC,FEDD,
	1              AQWEXT,HQWEXT,KQWEXT,NQWEXT,HMIDQWEXT,NMIDQWEXT,
	1              INBC,HBC_CMF(1),IPLUS,FL,dLOG_NU,DIF,DBB,IC,
	1              VDOP_VEC_EXT,DELV_FRAC_FG,REXT_FAC,
	1              METHOD,FG_SOL_OPTIONS,THK_CONT,
	1              FIRST_FREQ,NEW_FREQ,
	1              NCEXT,NPEXT,NDEXT)
	     END IF
!
	     CALL TUNE(ITWO,'FG_J_CMF_ACC')
!
! RSQHNU/TB not needed as computed by MOM routines.
!	     CALL UNGRID(RSQHNU,ND,TB,NDEXT,POS_IN_NEW_GRID)
!
! We assume Rayleigh scattering, when present, is coherent in the CMF.
! NB: We have to pass a non zero CHI_RAY_EXT to MOM_J_CMF_V6, since it is
! interpolated in the LOG plane.
!
	     IF(COHERENT_ES)THEN
	       TA(1:NDEXT)=ETAEXT(1:NDEXT)
	       TB(1:NDEXT)=CHI_SCAT_EXT(1:NDEXT)
	     ELSE
	       TA(1:NDEXT)=ETAEXT(1:NDEXT)+ESECEXT(1:NDEXT)*RJEXT_ES(1:NDEXT)
	       TB(1:NDEXT)=CHI_RAY_EXT(1:NDEXT)
	     END IF
	     CALL TUNE(IONE,'MOM_J_CMF_ACC')
             IF(PLANE_PARALLEL_NO_V)THEN
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling MOM_J_PP_V1 with accurate grid'
	       CALL MOM_J_PP_V1(TA,CHIEXT,TB,
	1                  REXT,FEDD,RJEXT,RSQHNU,HBC_CMF,NBC_CMF,INBC,
	1                  FL,DIF,DBB,IC,METHOD,L_TRUE,
	1                  IZERO,FIRST_FREQ,NEW_FREQ,NDEXT)
	       HFLUX_AT_OB=HBC_CMF(1)*RJ(1)-HBC_CMF(2)
	       IF(.NOT. DIF)HFLUX_AT_OB=0.5D0*IC*(0.5D0+INBC)-INBC*RJ(ND)
             ELSE IF(PLANE_PARALLEL)THEN
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling PP_MOM_CMF_V1 with accurate grid'
               CALL PP_MOM_CMF_V1(TA,CHIEXT,TB,VEXT,SIGMAEXT,REXT,
	1                  FEDD,GEDD,N_ON_J,RJEXT,RSQHNU,
	1                  VDOP_VEC_EXT,DELV_FRAC_MOM,
	1                  INBC,HBC_CMF(1),HBC_CMF(2),NBC_CMF(1),NBC_CMF(2),
	1                  FL,dLOG_NU,DIF,DBB,IC,
	1                  N_TYPE,METHOD,L_TRUE,
	1                  FIRST_FREQ,NEW_FREQ,NDEXT)
	       HFLUX_AT_OB=HBC_CMF(1)*RJ(1)-HBC_CMF(2)
	       IF(.NOT. DIF)HFLUX_AT_OB=0.5D0*IC*(0.5D0+INBC)-INBC*RJ(ND)
	     ELSE IF(USE_DJDT_RTE)THEN
	      IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling MOM_J_DDT_V4'
	      CALL MOM_J_DDT_V4(TA,CHIEXT,TB,
	1              VEXT,REXT,FEDD,RJEXT,RSQHNU,XM,                  !Use XM for DJDt_TERM,
	1              HFLUX_AT_IB,HFLUX_AT_OB,
	1              VDOP_VEC_EXT,DELV_FRAC_MOM,FL,dLOG_NU,DBB,
	1              H_CHK_OPTION,INNER_BND_METH,OUTER_BND_METH,
	1              METHOD,COHERENT_ES,FIRST_FREQ,NEW_FREQ,
	1              INCL_DJDT_TERMS,USE_DR4JDT,DJDT_RELAX_PARAM,NCEXT,NPEXT,NDEXT,NCF)

             ELSE IF(USE_J_REL)THEN
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling MOM_JREL_V8 with accurate grid'
	       CALL MOM_JREL_V8(TA,CHIEXT,TB,VEXT,SIGMAEXT,REXT,
	1             RJEXT,RSQHNU,HFLUX_AT_IB,HFLUX_AT_OB,
	1             VDOP_VEC_EXT,DELV_FRAC_MOM,FL,dLOG_NU,DBB,H_CHK_OPTION,
	1             IB_STAB_FACTOR,INNER_BND_METH,OUTER_BND_METH,
	1             METHOD,L_TRUE,N_TYPE,
	1             INCL_ADVEC_TERMS_IN_TRANS_EQ,INCL_REL_TERMS,FIRST_FREQ,NDEXT)
                 IF(LST_ITERATION)THEN
                   DO I=1,NDEXT
                     TA(I)=RJEXT(I)*REXT(I)*REXT(I)
                   END DO
                   T1=HFLUX_AT_IB*REXT(ND)*REXT(ND)
                   T2=HFLUX_AT_OB/RJEXT(1)
                   CALL OUT_JH(TA,RSQHNU,T1,T2,FL,NCF,REXT,VEXT,NDEXT,FIRST_FREQ,'NORMAL')
	        END IF
	     ELSE IF(USE_LAM_ES)THEN
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Using J from formal-solution for ray-by-ray solution'
	       RJEXT(1:NDEXT)=TC(1:NDEXT); RSQHNU(1:NDEXT)=0.0D0
	       CALL GET_RSQH_REL(RSQHNU, REXT, VEXT, FL, NDEXT)
	       IF(L .LT. 3)FOLD(1:ND)=0.0D0
	     ELSE
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling MOM_J_CMF_V8 with accurate grid'
	       CALL MOM_J_CMF_V8(TA,CHIEXT,TB,VEXT,SIGMAEXT,REXT,
	1  	       RJEXT,RSQHNU,VDOP_VEC_EXT,DELV_FRAC_MOM,
	1              FL,dLOG_NU,DIF,DBB,IC,
	1              N_TYPE,METHOD,L_TRUE,IONE,           !IONE => old outer bc
	1              FIRST_FREQ,NEW_FREQ,
	1              NCEXT,NPEXT,NDEXT)
	       IF(.NOT. DIF)HFLUX_AT_IB=0.5D0*IC*(0.5D0+INBC)-INBC*RJ(ND)
	       HFLUX_AT_OB=HBC_CMF(1)*RJ(1)
	     END IF
	     CALL TUNE(ITWO,'MOM_J_CMF_ACC')
!
! We set NEW_FREQ to false so that FG_J_CMF continues to use the same
! AV_PREV and CV_PREV. NEW_FREQ must be set to true again outside the
! F iteration loop.
!
	     NEW_FREQ=.FALSE.
!
! Update "inaccurate" iteration counter
!
	      L=L+1
!
! Check if F has converged.
!
	      INACCURATE=.FALSE.
	      IF(L .LT. 20 .OR. COMPUTE_EDDFAC)THEN
	        T1=0.0D0
	        DO I=1,NDEXT
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
!
	      IF(L .GT. 10)THEN
	         WRITE(LUER,*)'Possible error converging f - T1 is',T1
	         WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION
	      	 INACCURATE=.FALSE.
	      END IF
	    END DO
!
! Output RJ for subsequent iterations.
!
	    WRITE(LU_EDD,REC=ACCESS_F)(RJEXT(I),I=1,NDEXT),FL
!
! Store J on the normal mesh. No interpolation is involved here as
! it is assumed that the fine grid was created by the addition of extra
! points only.
!
	    CALL UNGRID(RJ,ND,RJEXT,NDEXT,POS_IN_NEW_GRID)
	    CALL UNGRID(RJ_ES,ND,RJEXT_ES,NDEXT,POS_IN_NEW_GRID)
	    CALL UNGRID(K_MOM,ND,FEDD,NDEXT,POS_IN_NEW_GRID)
!
! Compute K for use in computing mechanical energy loss.
!
	    K_MOM(1:ND)=RJ(1:ND)*K_MOM(1:ND)
!
	    IF(COHERENT_ES)THEN
	      SOURCE(1:ND)=ZETA(1:ND)+CHI_SCAT(1:ND)*RJ(1:ND)/CHI(1:ND)
	    ELSE
	      SOURCE(1:ND)=ZETA(1:ND) +
	1          (ESEC(1:ND)*RJ_ES(1:ND)+CHI_RAY(1:ND)*RJ(1:ND))/CHI(1:ND)
	    END IF
!
! Update record for next frequency
!
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
	    WRITE(LU_JCOMP,'(I7,ES14.6,2ES12.4,F10.2,3ES12.4,F10.2)')
	1                       FREQ_INDX,FL,
	1                       RJ(1),TC(1),T1,HBC_CMF(1),
	1                       RJ(ND),TC(NDEXT),T2
	    IF(ML .EQ. NCF)CLOSE(UNIT=LU_JCOMP)
	  END IF
	  IF(ML .EQ. NCF)THEN
	     WRITE(LUER,*)'Average number of calls to FG_J_CMF is',FG_COUNT/NCF
	  END IF
!
! Set up for the compuation of the observes flux. LST_ITERATION is
! TRUE if FLUX_CAL_ONLY is true (single iteration with coherent,
! last iteration if non-coherent).
!
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
	1         MAXCH .LT. VAL_DO_NG .AND. .NOT. USE_FORMAL_REL)  )THEN
!
! Quick and dirty method to ge an extended DENSITY vector. Will use TB in
! the call to CMF_FORM_SOL.
!
	       CALL EXTEND3OPAC(TA,TB,TC,COEF,INDX,NDEXT,
	1              DENSITY,DENSITY,DENSITY,ND)
!
	      IF(COHERENT_ES)THEN
	        TA(1:NDEXT)=ETAEXT(1:NDEXT)+ESECEXT(1:NDEXT)*RJEXT(1:NDEXT)
	      ELSE
	        TA(1:NDEXT)=ETAEXT(1:NDEXT)+ESECEXT(1:NDEXT)*RJEXT_ES(1:NDEXT)+
	1           CHI_RAY_EXT(1:NDEXT)*RJEXT(1:NDEXT)
	      END IF
!
! NB Using TA for ETA, U for P_OBS (temporay measure), I for NP_OBS.
!
	      CALL TUNE(IONE,'CMF_FORM_SOL')
	      CALL CMF_FORM_SOL(TA,CHIEXT,CHI_SCAT_EXT,
	1                 TB,VEXT,SIGMAEXT,REXT,PEXT,
	1                 P_OBS,IPLUS,NP_OBS,NP_OBS_MAX,
	1                 MU_AT_RMAX,HQW_AT_RMAX,RMAX_OBS,V_AT_RMAX,
	1                 FL,dLOG_NU,DIF,DBB,IC,METHOD,
	1                 EXTEND_FRM_SOL,INSERT_FREQ_FRM_SOL,
	1                 FRAC_DOP,V_DOP,dV_CMF_PROF,dV_CMF_WING,
	1                 FIRST_FREQ,NCEXT,NPEXT,NDEXT)
	      CALL TUNE(ITWO,'CMF_FORM_SOL')
	    ELSE IF(FIRST_FREQ)THEN
!
! So as defined for normal OBSFLUX calculation.
!
	      NP_OBS=NPEXT
	      P_OBS(1:NPEXT)=PEXT(1:NPEXT)
	      RMAX_OBS=R(1)
	      V_AT_RMAX=V(1)
	    END IF
!
	  CALL TUNE(ITWO,'CONT_VEL')
!
! 
!
	ELSE IF(CONT_VEL)THEN
!
! NB: CHI_PREV is used to refer to the continuum opacity at the previous
! frequency. Is does not need to be multiplied by CLUMP_FAC, as it is compared
! directly to CHI_CONT.
!
! For HBC and NBC only the first vector element is used.
!
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
!
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
	    IF(ABS(T1-FL)/FL .GT. 1.0D-14)THEN
	      WRITE(LUER,*)'Error - 4 incorrect reading of'//
	1                  ' the mean intensity'
	      WRITE(LUER,*)'Intensity is ',FL,'Old Intensity is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP
	    END IF
	  END IF
!
! If we are using incoherent electron scattering, RJEXT_ES must be available.
!
	  IF(.NOT. COHERENT_ES)THEN
	    READ(LU_ES,REC=ACCESS_F)(RJ_ES(I),I=1,ND),T1
	    IF(ABS(T1-FL)/FL .GT. 1.0D-14)THEN
	      WRITE(LUER,*)'Error - 5 incorrect reading of'//
	1                ' the mean intensity'
	      WRITE(LUER,*)'Intensity is ',FL,'Old Intensity is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP
	    END IF
	  END IF
!
! We will do this twice, so that F is of higher accuracy.
!
	  INACCURATE=.TRUE.
	  L=0
	  DO WHILE(INACCURATE)
!
	     IF(COHERENT_ES)THEN
	       TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ(1:ND)
	     ELSE
	       TA(1:ND)=ETA_CLUMP(1:ND)+ESEC_CLUMP(1:ND)*RJ_ES(1:ND)+
	1                                CHI_RAY_CLUMP(1:ND)*RJ(1:ND)
	     END IF
!
! NB Using TA for ETA, TC for JNU_VEC, and TB for HNU_VEC
!
	     CALL TUNE(IONE,'FG_J_CMF')
	     FG_COUNT=FG_COUNT+1
	     IF(PLANE_PARALLEL_NO_V)THEN
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling FCOMP_PP'
	       SOURCE(1:ND)=TA(1:ND)/CHI_CLUMP(1:ND)
	       CALL FCOMP_PP_V2(R,TC,FEDD,SOURCE,CHI_CLUMP,IPLUS,HBC_CMF,
	1               NBC_CMF,INBC,DBB,IC,THK_CONT,DIF,ND,NC,METHOD)
	     ELSE IF(PLANE_PARALLEL)THEN
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling PP_FORM_CMF_V2'
	       CALL PP_FORM_CMF_V2(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,
	1               TC,TB,FEDD,GEDD,N_ON_J,INBC,
	1               HBC_CMF(1),HBC_CMF(2),NBC_CMF(1),NBC_CMF(2),
	1               IPLUS,FL,dLOG_NU,DIF,DBB,IC,VDOP_VEC,DELV_FRAC_FG,
	1               METHOD,FG_SOL_OPTIONS,THK_CONT,L_FALSE,
	1               FIRST_FREQ,NEW_FREQ,N_TYPE,NC,ND)
             ELSE IF(USE_FORMAL_REL)THEN
               IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling CMF_FORMAL_REL_V4 in COMP_JCONT_V4.INC'
               CALL CMF_FORMAL_REL_V4
	1                 (TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,P,
	1                  TC,FEDD,HFLUX_AT_IB,HFLUX_AT_OB,IPLUS,
	1                  FL,dLOG_NU,BNUE,DBB,
	1                  INNER_BND_METH,THK_CONT,
	1                  VDOP_VEC,DELV_FRAC_FG,REXT_FAC,
	1                  METHOD,FIRST_FREQ,NEW_FREQ,NC,NP,ND)
	       INBC=0.0D0
	     ELSE
	        IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling FG_J_CMF_V12'
	        CALL FG_J_CMF_V12(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,P,
	1                  TC,FEDD,AQW,HQW,KQW,NQW,HMIDQW,NMIDQW,
	1                  INBC,HBC_CMF(1),IPLUS,FL,dLOG_NU,DIF,DBB,IC,
	1                  VDOP_VEC,DELV_FRAC_FG,REXT_FAC,
	1                  METHOD,FG_SOL_OPTIONS,THK_CONT,
	1                  FIRST_FREQ,NEW_FREQ,NC,NP,ND)
	     END IF
!
	     CALL TUNE(ITWO,'FG_J_CMF')
!
! We assume that the Rayleig scattering component is coherent.
! Only the e.s. component may be incoherent. NB: We have to pass
! a non zero CHI_RAY to MOM_J_CMF_V6, since it is interpolated
! in the LOG plane.
!
	     IF(COHERENT_ES)THEN
	       TA(1:ND)=ETA_CLUMP(1:ND)
	       TB(1:ND)=CHI_SCAT_CLUMP(1:ND)
	     ELSE
	       TA(1:ND)=ETA_CLUMP(1:ND)+ESEC_CLUMP(1:ND)*RJ_ES(1:ND)
               TB(1:ND)=CHI_RAY_CLUMP(1:ND)
	     END IF
	     CALL TUNE(IONE,'MOM_J_CMF')
             IF(PLANE_PARALLEL_NO_V)THEN
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling MOM_J_PP_V1'
	       CALL MOM_J_PP_V1(TA,CHI_CLUMP,TB,
	1                  R,FEDD,RJ,RSQHNU,HBC_CMF,NBC_CMF,INBC,
	1                  FL,DIF,DBB,IC,METHOD,L_TRUE,
	1                  IZERO,FIRST_FREQ,NEW_FREQ,ND)
	       HFLUX_AT_OB=HBC_CMF(1)*RJ(1)-HBC_CMF(2)
	       IF(.NOT. DIF)HFLUX_AT_OB=0.5D0*IC*(0.5D0+INBC)-INBC*RJ(ND)
             ELSE IF(PLANE_PARALLEL)THEN
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling PP_MOM_CMF_V1'
	       CALL PP_MOM_CMF_V1(TA,CHI_CLUMP,TB,V,SIGMA,R,
	1                  FEDD,GEDD,N_ON_J,RJ,RSQHNU,
	1                  VDOP_VEC,DELV_FRAC_MOM,
	1                  INBC,HBC_CMF(1),HBC_CMF(2),NBC_CMF(1),NBC_CMF(2),
	1                  FL,dLOG_NU,DIF,DBB,IC,
	1                  N_TYPE,METHOD,L_TRUE,
	1                  FIRST_FREQ,NEW_FREQ,ND)
               HFLUX_AT_OB=HBC_CMF(1)*RJ(1)-HBC_CMF(2)
               IF(.NOT. DIF)HFLUX_AT_OB=0.5D0*IC*(0.5D0+INBC)-INBC*RJ(ND)
!
	    ELSE IF(USE_DJDT_RTE)THEN
	      IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling MOM_J_DDT_V3'
	      CALL MOM_J_DDT_V4(TA,CHI_CLUMP,TB,
	1              V,R,FEDD,RJ,RSQHNU,XM,                    !Use TC for DJDt_TERM,
	1              HFLUX_AT_IB,HFLUX_AT_OB,
	1              VDOP_VEC,DELV_FRAC_MOM,FL,dLOG_NU,DBB,
	1              H_CHK_OPTION,INNER_BND_METH,OUTER_BND_METH,
	1              METHOD,COHERENT_ES,FIRST_FREQ,NEW_FREQ,
	1              INCL_DJDT_TERMS,USE_DR4JDT,DJDT_RELAX_PARAM,NC,NP,ND,NCF)

             ELSE IF(USE_J_REL)THEN
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling MOM_JREL_V7'
	       CALL MOM_JREL_V8(TA,CHI_CLUMP,TB,V,SIGMA,R,
	1             RJ,RSQHNU,HFLUX_AT_IB,HFLUX_AT_OB,VDOP_VEC,DELV_FRAC_MOM,
	1             FL,dLOG_NU,DBB,H_CHK_OPTION,
	1             IB_STAB_FACTOR,INNER_BND_METH,OUTER_BND_METH,
	1             METHOD,L_TRUE,N_TYPE,
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
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Using J from formal-solution for ray-by-ray solution'
	       RJ(1:ND)=TC(1:ND); RSQHNU(1:ND)=0.0D0
	       CALL GET_RSQH_REL(RSQHNU, R, V, FL, ND)
	       IF(L .LT. 3)FOLD(1:ND)=0.0D0
	     ELSE
	       IF(FG_COUNT .EQ. 1)WRITE(LUER,*)'Calling MOM_J_CMF_V8'
	       CALL MOM_J_CMF_V8(TA,CHI_CLUMP,TB,V,SIGMA,R,
	1  	       RJ,RSQHNU,VDOP_VEC,DELV_FRAC_MOM,
	1              FL,dLOG_NU,DIF,DBB,IC,
	1              N_TYPE,METHOD,L_TRUE,IONE,
	1              FIRST_FREQ,NEW_FREQ,NC,NP,ND)
               IF(.NOT. DIF)HFLUX_AT_OB=0.5D0*IC*(0.5D0+INBC)-INBC*RJ(ND)
               HFLUX_AT_OB=HBC_CMF(1)*RJ(1)
	     END IF
	     CALL TUNE(ITWO,'MOM_J_CMF')
!
! We set NEW_FREQ to false so that FG_J_CMF continues to use the same
! AV_PREV and CV_PREV. NEW_FREQ must be set to true again outside the
! F iteration loop.
!
	     NEW_FREQ=.FALSE.
!
! Update "inaccurate" iteration counter
!
	      L=L+1
!
! Check if F has converged.
!
	      INACCURATE=.FALSE.
	      IF(L .LT. 20 .OR. COMPUTE_EDDFAC)THEN
	        T1=0.0D0
	        DO I=1,ND
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
!
	      IF(L .GT. 10)THEN
	         WRITE(LUER,*)'Possible error converging f - T1 is',T1
	         WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION
	      	 INACCURATE=.FALSE.
	      END IF
	    END DO
!
	    IF(COHERENT_ES)THEN
	      SOURCE(1:ND)=ZETA(1:ND)+CHI_SCAT(1:ND)*RJ(1:ND)/CHI(1:ND)
	    ELSE
	      SOURCE(1:ND)=ZETA(1:ND) +
	1          (ESEC(1:ND)*RJ_ES(1:ND)+CHI_RAY(1:ND)*RJ(1:ND))/CHI(1:ND)
	    END IF
!
! Output RJ for subsequent iterations.
!
	    WRITE(LU_EDD,REC=ACCESS_F)(RJ(I),I=1,ND),FL
!
! Compute K for use in computing mechanical energy loss.
!
	    DO I=1,ND
	      K_MOM(I)=RJ(I)*FEDD(I)
	    END DO
!
! Update record for next frequency
!
	    ACCESS_F=ACCESS_F+1
	    FL_OLD=FL
!
! Set up for the compuation of the observes flux. LST_ITERATION is
! TRUE if FLUX_CAL_ONLY is true (single iteration with coherent,
! last iteration if non-coherent).
!
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
	1         MAXCH .LT. VAL_DO_NG .AND. .NOT. USE_FORMAL_REL) )THEN
	      IF(COHERENT_ES)THEN
     	        TA(1:ND)=ETA_CLUMP(1:ND)+CHI_SCAT_CLUMP(1:ND)*RJ(1:ND)
	      ELSE
     	        TA(1:ND)=ETA_CLUMP(1:ND)+ESEC_CLUMP(1:ND)*RJ_ES(1:ND) +
	1                                  CHI_RAY_CLUMP(1:ND)*RJ(1:ND)
	      END IF
!
! NB Using TA for ETA, U for P_OBS (temporay measure), I for NP_OBS.
!
	      CALL TUNE(IONE,'CMF_FORM_SOL')
	      CALL CMF_FORM_SOL(TA,CHI_CLUMP,CHI_SCAT_CLUMP,
	1                 DENSITY,V,SIGMA,R,P,
	1                 P_OBS,IPLUS,NP_OBS,NP_OBS_MAX,
	1                 MU_AT_RMAX,HQW_AT_RMAX,RMAX_OBS,V_AT_RMAX,
	1                 FL,dLOG_NU,DIF,DBB,IC,METHOD,
	1                 EXTEND_FRM_SOL,INSERT_FREQ_FRM_SOL,
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
!
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
	    WRITE(LU_JCOMP,'(I7,ES14.6,2ES12.4,F10.2,3ES12.4,F10.2)')
	1                       FREQ_INDX,FL,
	1                       RJ(1),TC(1),T1,HBC_CMF(1),
	1                       RJ(ND),TC(ND),T2
	    IF(ML .EQ. NCF)CLOSE(UNIT=LU_JCOMP)
	  END IF
	  IF(ML .EQ. NCF)THEN
	     WRITE(LUER,*)'Average number of calls to FG_J_CMF is',FG_COUNT/NCF
	  END IF
	  CALL TUNE(ITWO,'CONT_VEL')
!
! 
!
	ELSE IF(EDDINGTON)THEN
!
! Calculation of "static" J in the continuum using Edington factors.
!
	  CALL TUNE(IONE,'JFEAU')
	  IF(COMPUTE_EDDFAC)THEN
	    DO I=1,ND
	      RJ(I)=0.0D0
              FOLD(I)=FEDD(I)
	    END DO
	  ELSE
	    READ(LU_EDD,REC=ACCESS_F)(RJ(I),I=1,ND),T1
	    IF(ABS(T1-FL)/FL .GT. 1.0D-14)THEN
	      WRITE(LUER,*)'Error - incorrect reading of Mean Intensity'
	      WRITE(LUER,*)'Intensity is ',FL,'Old Intensity is ',T1
	      WRITE(LUER,*)'Error occurred in '//SECTION
	      STOP
	    END IF
	  END IF
!
! We will do this twice, so that F is of higher accuracy.
!
	  INACCURATE=.TRUE.
	  L=0
	  DO WHILE(INACCURATE)
	    DO I=1,ND
	      SOURCE(I)=ZETA(I)+THETA(I)*RJ(I)
	    END DO
	    S1=SOURCE(1)
!
	    CALL FQCOMP_IBC_V2(TA,TB,TC,XM,DTAU,R,Z,P,QEDD,FEDD,
	1            SOURCE,CHI_CLUMP,dCHIdR,AQW,KQW,
	1            DBB,HBC_J,HBC_S,INBC,IC,
	1            THK_CONT,INNER_BND_METH,NC,ND,NP,METHOD)
	    CALL JFEAU_IBC_V2(TA,TB,TC,DTAU,R,RJ,QEDD,FEDD,
	1          ZETA,THETA,CHI_CLUMP,DBB,IC,HBC_J,HBC_S,
	1          INBC,THK_CONT,INNER_BND_METH,ND,METHOD)
!
! Update "inaccurate" iteration counter
!
	      L=L+1
!
! Check if F has converged.
!
	      INACCURATE=.FALSE.
	      IF(L .LT. 3 .OR. COMPUTE_EDDFAC)THEN
	        T1=0.0D0
	        DO I=1,ND
	          T1=MAX(ABS(FOLD(I)-FEDD(I)),T1)
	          FOLD(I)=FEDD(I)
	        END DO
	        IF(T1 .GT. ACC_EDD_FAC)INACCURATE=.TRUE.
	      END IF
!
	      L=L+1
	      IF(L .GT. 15)THEN
	         WRITE(LUER,*)'Possible error converging f - T1 is',T1
	         WRITE(LUER,*)'Frequency is ',FL,' in section '//SECTION
	      	 INACCURATE=.FALSE.
	      END IF

	    END DO
!
	    DO I=1,ND
	      K_MOM(I)=RJ(I)*FEDD(I)
	    END DO
!
! Output mean intensity for subsequent iterations.
!
	    WRITE(LU_EDD,REC=ACCESS_F)(RJ(I),I=1,ND),FL
!
! Update record for next frequency
	    ACCESS_F=ACCESS_F+1
!
	  CALL TUNE(ITWO,'JFEAU')
	ELSE
!
! Calculation of "static" J in the continuum using Rybick method.
!
	  CALL TUNE(IONE,'JSOL')
	  CALL NEWJSOLD(TA,TB,TC,XM,WM,FB,RJ,DTAU,R,Z,P,
	1       ZETA,THETA,CHI_CLUMP,dCHIdR,AQW,
	1       THK_CONT,DIF,DBB,IC,NC,ND,NP,METHOD)
	  CALL TUNE(ITWO,'JSOL')
!
! Compute K_MOM.
!
	  DO I=1,ND
	    SOURCE(I)=ZETA(I)+THETA(I)*RJ(I)
	  END DO
	  S1=SOURCE(1)
	  CALL FQCOMP_IBC(TA,TB,TC,XM,DTAU,R,Z,P,QEDD,FEDD,
	1            SOURCE,CHI_CLUMP,dCHIdR,AQW,KQW,DBB,HBC_J,HBC_S,
	1            INBC,IC,THK_CONT,DIF,NC,ND,NP,METHOD)
	  DO I=1,ND
	    K_MOM(I)=RJ(I)*FEDD(I)
	  END DO
	END IF
!
	IF(ABS(RJ(1)) .GT. 1.0D+30)THEN
	  WRITE(LUER,*)' '
	  WRITE(LUER,*)' '
	  WRITE(LUER,*)('*',I=1,80)
	  WRITE(LUER,*)'Error in COMP_JCONT_V4.INC'
	  WRITE(LUER,*)'Mean intensity blowing up, which is due to an instabilty'
	  WRITE(LUER,*)'Try a diferent N_TYPE option'
	  WRITE(LUER,*)'Alternatively set INC_GRID=ALL_FREQ=T and NINS=1'
	  WRITE(LUER,*)'Setting DELV_FG & DELV_MOM to 1 will also help'
	  WRITE(LUER,*)('*',I=1,80)
	  WRITE(LUER,*)' '
	  WRITE(LUER,*)' '
	  STOP
	END IF
