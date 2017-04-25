	SUBROUTINE PLTPHOT_SUB(XSPEC,XV,YV,WV,ZV,NF_MAX,CROSS,NMAX,XAXIS,YAXIS,XRAYS,VSM_DIE_KMS,ND)
	USE MOD_DISP
	USE MOD_USR_OPTION
	USE MOD_USR_HIDDEN
	USE MOD_WR_STRING
	USE MOD_LEV_DIS_BLK
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered: 19-Jun-2016 - Bug fix. Routine was returning to zero cross section for NU > NU(edge)
!                           when ED_VAL was non zero (i.e., level dissolution switched on).
!                           Level dissolution only used when plotting individual cross-sections.
	INTEGER NMAX
	INTEGER NF_MAX
	INTEGER ND
	REAL*8 XV(NF_MAX)
	REAL*8 YV(NMAX)
	REAL*8 WV(NMAX)
	REAL*8 ZV(NMAX)
	REAL*8 CROSS(NMAX)
	REAL*8 VSM_DIE_KMS
!
	CHARACTER(LEN=*) XAXIS
	CHARACTER(LEN=*) YAXIS
	LOGICAL XRAYS
!
	REAL*8 FREQ
	REAL*8 FREQ_MAX
	REAL*8 FREQ_RES
	REAL*8 EDGE_FREQ
	REAL*8 ED_VAL
	REAL*8 MIN_ED,MAX_ED
	REAL*8 T1,T2,T3
	REAL*8 EXC_EN
	REAL*8 DIS_CONST
	REAL*8 XDIS,YDIS
        REAL*8 EV_TO_HZ,ANG_TO_HZ
	REAL*8 C_CMS
	REAL*8 TMP_GION
!
	INTEGER I,J,K
	INTEGER CNT
	INTEGER ID
	INTEGER ISPEC
	INTEGER PHOT_ID
	LOGICAL DO_SPECIES
	LOGICAL FLAG
	LOGICAL DIE_REG
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_OUT=6
	LOGICAL,PARAMETER :: L_FALSE=.FALSE.
	LOGICAL DO_EV
!
	CHARACTER(LEN=10) XSPEC
        CHARACTER(LEN=120) DEFAULT
!
	REAL*8 XCROSS_V2
	REAL*8 SPEED_OF_LIGHT
	CHARACTER*30 UC
	EXTERNAL SUB_PHOT_GEN
	EXTERNAL UC,SPEED_OF_LIGHT
	EXTERNAL XCROSS_v2
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Conversion factor from Kev to units of 10^15 Hz.
! Conversion factor from Angstroms to units of 10^15 Hz.
!
	DO_EV=.TRUE.
        EV_TO_HZ=0.241838D+00
        ANG_TO_HZ=SPEED_OF_LIGHT()*1.0D-07      !10^8/10^15
        C_CMS=SPEED_OF_LIGHT()
!
	WRITE(6,*)NSPEC,XSPEC
	DO_SPECIES=.FALSE.
	DO ISPEC=1,NSPEC
	  IF(XSPEC .EQ. UC(SPECIES(ISPEC)))THEN
	    DO_SPECIES=.TRUE.
	    IF(SPECIES_BEG_ID(ISPEC) .LE.  0)THEN
	      WRITE(T_OUT,*)'Error --- this  species is unavailable'
	      RETURN
	    END IF
	    EXIT
	  END IF
	END DO
!
	CALL USR_OPTION(PHOT_ID,'PHOT_ID','1','Photoionization route')
	FREQ_RES=MIN(3000.0D0,VSM_DIE_KMS)/2.0D0
	DEFAULT=WR_STRING(FREQ_RES)
	CALL USR_OPTION(FREQ_RES,'FREQ_RES',DEFAULT,'Frequency resolution in km/s')
	FREQ_RES=FREQ_RES/3.0D+05
!
	IF(DO_SPECIES)THEN
	  XAXIS='\gn(10\u15 \dHz)'
	  FREQ_MAX=1000.0D0
	  DEFAULT=WR_STRING(FREQ_MAX)
	  CALL USR_OPTION(FREQ_MAX,'FREQ_MAX',DEFAULT,'Maximum frequency in units of 10^15 Hz')
	  DEFAULT='T'
	  CALL USR_OPTION(DO_EV,'DO_EV',DEFAULT,'Plot X-axis in eV?')
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    FREQ=ATM(SPECIES_BEG_ID(ISPEC))%EDGEXzV_F(ID)
	    J=0
	    DO WHILE(FREQ .LT. FREQ_MAX)
	      IF(J .EQ. NF_MAX)EXIT
	      J=J+1
	      ZV(J)=FREQ
	      IF(FREQ .GE. ATM(ID)%EDGEXzV_F(1))THEN
	        CALL SUB_PHOT_GEN(ID,CROSS,FREQ,ATM(ID)%EDGEXzV_F,ATM(ID)%NXzV_F,PHOT_ID,FLAG)
	        YV(J)=CROSS(1)
	      ELSE
	        YV(J)=0.0D0
	      END IF
	      IF(XRAYS .AND. ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
                T2=AT_NO(SPECIES_LNK(ID))+1-ATM(ID)%ZXzV
                T1=XCROSS_V2(FREQ,AT_NO(SPECIES_LNK(ID)),T2,IZERO,IZERO,L_FALSE,L_FALSE)
	        YV(J)=YV(J)+T1
	      END IF
	      YV(J)=1.0D+08*YV(J)
	      FREQ=FREQ*(1.0D0+FREQ_RES)
	    END DO
	    IF(DO_EV)THEN
	      ZV(1:J)=ZV(1:J)/EV_TO_HZ
	      XAXIS='\gn(eV)'
	    END IF
	    CALL DP_CURVE(J,ZV,YV)
	    WRITE(6,'(A,I2,A,A)')' Curve ',ID-SPECIES_BEG_ID(ISPEC)+1,' is: ',TRIM(ION_ID(ID))
	  END DO
	  YAXIS='\gs(Mb)'
	  RETURN
	END IF
!
	CALL USR_OPTION(I,'LEV','1','Level ID: Use WRID to check levs')
	ED_VAL=0.0D0
	WRITE(T_OUT,'(A)')' '
	WRITE(T_OUT,'(A)')' Input non-zero Ne for level-dissolution.'
	WRITE(T_OUT,'(A)')' ED is resticted to the range covered by the model atmosphere.'
	WRITE(T_OUT,'(A)')' Closest model atmospere value is used'
	WRITE(T_OUT,'(A)')' '
	IF(PHOT_ID .EQ. 1)CALL USR_OPTION(ED_VAL,'ED_VAL','0.0D0','Approximate Ne value')
	IF(ED_VAL .GT. 0)THEN
	  MIN_ED=MINVAL(ED)
	  MAX_ED=MAXVAL(ED)	
	  IF(ED_VAL .LT. MIN_ED)THEN
	    ED_VAL=MIN_ED
	    CNT=MINLOC(ED,IONE)
	  ELSE IF(ED_VAL .GT. MAX_ED)THEN
	    ED_VAL=MAX_ED
	    CNT=MAXLOC(ED,IONE)
	  ELSE
	    DO J=1,ND-1
	     IF(ED_VAL .GT. ED(J) .AND. ED_VAL .LE. ED(J+1))THEN
	       CNT=J
	       IF( LOG(ED_VAL/ED(J)) .GT. LOG(ED(J+1)/ED_VAL))CNT=J+1
	       EXIT
	     END IF
	    END DO
	  END IF
	END IF
!
	FREQ_RES=MIN(3000.0D0,VSM_DIE_KMS)/2.0D0
	DEFAULT=WR_STRING(FREQ_RES)
	CALL USR_HIDDEN(FREQ_RES,'FREQ_RES',DEFAULT,'Frequency resolution in km/s')
	FREQ_RES=FREQ_RES/3.0D+05
!
	DO ID=1,NUM_IONS
	  IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	    IF(I .GT. ATM(ID)%NXzV_F)THEN
	      WRITE(T_OUT,*)'Invalid level ID for this species'
	      RETURN
	    ELSE
	      WRITE(T_OUT,*)RED_PEN
	      WRITE(T_OUT,*)'Level is: ',ATM(ID)%XzVLEVNAME_F(I)
	      WRITE(T_OUT,'(A,F11.8)')'Ionization energy to g.s. (in 10^15 Hz) is: ',ATM(ID)%EDGEXzV_F(I)
	      IF(ED_VAL .NE. 0.0D0)THEN
	        WRITE(T_OUT,'(A,ES10.2)')' Photoionization cross-section evaluated at Ne=',ED(CNT)
	      END IF
	      WRITE(T_OUT,*)DEF_PEN
	    END IF
	    FLAG=.FALSE.				!Don't return edge value.
	    IF(ED_VAL .NE. 0.0D0)FLAG=.TRUE.
	    FREQ=ATM(ID)%EDGEXzV_F(I)
	    IF(FLAG)FREQ=0.7D0*ATM(ID)%EDGEXzV_F(I)
	    J=0
	    FREQ_MAX=20.0D0*ATM(ID)%EDGEXzV_F(I)
	    DEFAULT=WR_STRING(FREQ_MAX)
	    CALL USR_HIDDEN(FREQ_MAX,'FREQ_MAX',DEFAULT,'Maximum frequency in units of 10^15 Hz')
	    DO WHILE(FREQ .LT. FREQ_MAX)
	      CALL SUB_PHOT_GEN(ID,CROSS,FREQ,ATM(ID)%EDGEXzV_F,ATM(ID)%NXzV_F,PHOT_ID,FLAG)
	      IF(FLAG .AND. FREQ .LT. ATM(ID)%EDGEXzV_F(I))THEN
	        T1=ATM(ID)%ZXzV**3
	        T2=SQRT(3.289395*ATM(ID)%ZXzV*ATM(ID)%ZXzV/(ATM(ID)%EDGEXzV_F(I)-FREQ))
	        IF(T2 .GT. 2*ATM(ID)%ZXzV)THEN
	          T3=MIN(1.0D0,16.0D0*T1/(1.0D0+T2)/(1.0D0+T2)/3.0D0)
	          DIS_CONST=( T3*ATM(ID)%ZXzV*T1/(T2**4) )**1.5D0
	          K=CNT
	          YDIS=1.091*(X_LEV_DIS(K)+4.0D0*(ATM(ID)%ZXzV-1)*A_LEV_DIS(K))*
	1                         B_LEV_DIS(K)*B_LEV_DIS(K)
	          XDIS=B_LEV_DIS(K)*X_LEV_DIS(K)
	          T1=7.782+XDIS*DIS_CONST
	          T2=T1/(T1+YDIS*DIS_CONST*DIS_CONST)
	          CROSS(I)=CROSS(I)*T2
	        ELSE
	          CROSS(I)=0.0D0
	        END IF
	      END IF
!
	      IF(XRAYS .AND. ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
                T2=AT_NO(SPECIES_LNK(ID))+1-ATM(ID)%ZXzV
                T1=XCROSS_V2(FREQ,AT_NO(SPECIES_LNK(ID)),T2,IZERO,IZERO,L_FALSE,L_FALSE)
	        CROSS(I)=CROSS(I)+T1
	      END IF
!
	      J=J+1
	      ZV(J)=FREQ/ATM(ID)%EDGEXzV_F(I)
	      YV(J)=1.0D+08*CROSS(I)
	      IF(J .EQ. NF_MAX)EXIT  
	      IF(FREQ .LT. ATM(ID)%EDGEXzV_F(I))THEN
	        FREQ=FREQ*(1.0D0+10.0D0/3.0D+05)
	      ELSE
	         FREQ=FREQ*(1.0D0+FREQ_RES)
	      END IF
	    END DO
	    EDGE_FREQ=ATM(ID)%EDGEXzV_F(I)
	    EXIT
	  END IF
	END DO
!

	CALL USR_OPTION(DIE_REG,'CUM','F','Plot recombination cummulative function?')
	IF(DIE_REG)THEN
	  CALL USR_OPTION(FREQ,'T','1.0','Input T (in 10^4 K)')
	  EXC_EN=0.0D0
          IF(PHOT_ID .NE. 1)THEN
	    CALL USR_OPTION(EXC_EN,'EXC_EN',' ','Excitaiton Energy (cm^-1) of final state')
	  END IF
	  EXC_EN=1.0D-15*C_CMS*EXC_EN
	  CALL USR_OPTION(TMP_GION,'GION',' ','G for ION (No def)')
!
	  T1=HDKT*ATM(ID)%EDGEXzV_F(I)/FREQ
	  WV(1:J)=0.0D0
	  T3=YV(1)*ZV(1)*ZV(1)*EXP(-T1*(ZV(1)-1.0D0))
	  DO K=2,J
	    T2=T3
	    T3=YV(K)*ZV(K)*ZV(K)*EXP(-T1*(ZV(K)-1.0D0))
	    WV(K)=WV(K-1)+0.5D0*(T2+T3)*(ZV(K)-ZV(K-1))
	  END DO
!
! Not that YV above is in Mbarns.
!
          T2=5.7885E-15*WV(J)*(ATM(ID)%EDGEXzV_F(I)**3)*ATM(ID)%GXzV_F(I)/TMP_GION
          T2=T2*EXP(-HDKT*EXC_EN/T1)/(FREQ**1.5)
          WRITE(6,*)'The recomination rate is:',T2
	  DO K=1,J
	    WV(K)=WV(K)/WV(J)
	  END DO  
	  XAXIS='\gn/\gn\do\u'
	  YAXIS='F(\gv)'
          YV(1:J)=WV(1:J)
	ELSE
	  XAXIS='\gn/\gn\do\u'
	  YAXIS='\gs(Mb)'
	END IF
!
	CALL USR_OPTION(FLAG,'LAM','F','Plot against wavelength?')
	IF(FLAG)THEN
	  DO K=1,J
	    ZV(K)=ANG_TO_HZ/(ZV(K)*ATM(ID)%EDGEXzV_F(I))
	  END DO
	  XAXIS='\gl(\A)'
	ELSE
	  CALL USR_OPTION(FLAG,'NNF','F','Plot against non-normalized frequency?')
	  IF(FLAG)THEN
	    ZV(1:J)=ZV(1:J)*EDGE_FREQ
	    XAXIS='\gn(10\u15 \dHz)'
	  END IF
	END IF
	CALL DP_CURVE(J,ZV,YV)
!
	RETURN
	END
