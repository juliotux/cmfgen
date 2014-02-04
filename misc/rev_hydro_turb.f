!
! Auxilary program designed to modify the HYDRO file output from CMFGEN.
! Progam modifies the adopted stellar mass. Not that the percentage error
! os now defined so that it has a rang of pm 200%.
!
	PROGRAM REV_HYDRO_FILE
	USE GEN_IN_INTERFACE
!
! Altered: 23-Jun-2003 - ND and STRING now initialized.
! Cleaned: 07-Nov-2000
!
	IMPLICIT NONE
!
	INTEGER I,J,IOS
	INTEGER NSTR
	INTEGER ND
	INTEGER, PARAMETER :: LU_OUT=11
	INTEGER, PARAMETER :: T_OUT=6
!
	REAL*8 R
	REAL*8 V
	REAL*8 E
	REAL*8 VdVdR
	REAL*8 dPdR
	REAL*8 g_tot
	REAL*8 g_rad
	REAL*8 g_elec
	REAL*8 Gamma
	REAL*8 MT
	REAL*8 RND
	REAL*8 MASS_NEW
	REAL*8 MASS_OLD
	REAL*8 GSUR_NEW
	REAL*8 GSUR_OLD
	REAL*8 GRAV_CON
!
	REAL*8 dTPdR
        REAL*8 VTURB
!
! These variables are used to estimate a new mass for the STAR so that the
! Hydrostatic equation is better satified (in a least squares sense) over
! some range of depths.
!
        REAL*8 SUM_ERROR,SUM_R,DENOM,DEL_M
        INTEGER LOW_LIM,HIGH_LIM
!
	REAL*8 GRAVITATIONAL_CONSTANT, MASS_SUN
	EXTERNAL GRAVITATIONAL_CONSTANT, MASS_SUN
!
	CHARACTER*132 STRING(500)
	CHARACTER*132 FMT
	CHARACTER*132 FILENAME
!
	REAL*8 P_VEL(200)
	REAL*8 P_REQ(200)
	REAL*8 P_GRAD(200)
	REAL*8 P_GELEC(200)
!
	GRAV_CON=1.0D-20*GRAVITATIONAL_CONSTANT()*MASS_SUN()
!
	ND=0
	STRING(:)=' '
	FILENAME='HYDRO'
	CALL GEN_IN(FILENAME,'Input hydro file')
	OPEN(UNIT=10,FILE=FILENAME,ACTION='READ',STATUS='OLD')
	  NSTR=0        
	  DO WHILE(1 .EQ. 1)
	    READ(10,'(A)',END=1000)STRING(NSTR+1)
	    IF(STRING(NSTR+1) .EQ. ' ' .AND. ND .EQ. 0)ND=NSTR-1
	    NSTR=NSTR+1
	  END DO         
	CLOSE(UNIT=10)
1000	CONTINUE
!
	WRITE(T_OUT,*)'Number of depth points is',ND
        I=160
!
	FILENAME='REV_HYDRO'
	CALL GEN_IN(FILENAME,'Output hydro file')
	CALL GEN_ASCI_OPEN(LU_OUT,FILENAME,'UNKNOWN',' ',' ',I,IOS)
	WRITE(LU_OUT,'(5X,A,15X,A,6X,A, 6(4X,A,1X), 6X,A, 6X,A)')
	1        'R','V','% Error','    VdVdR',
	1                   ' dPdR/ROH',
	1                   'dTPdR/ROH',
	1                   '    g_TOT',
	1                   '    g_RAD',
	1                    '   g_ELEC','Gamma','M(t)'
!
	READ(STRING(ND+1),*)RND
	DO I=1,NSTR
	  IF(INDEX(STRING(I),'Surface gravity is:') .NE. 0)THEN
	    J=INDEX(STRING(I),':')
	    READ(STRING(I)(J+1:),*)GSUR_OLD
	    MASS_OLD=GSUR_OLD/GRAV_CON*RND*RND
	    WRITE(T_OUT,'(1X,A,F8.2)')'Old mass is ',MASS_OLD
	    MASS_NEW=MASS_OLD
	    CALL GEN_IN(MASS_NEW,'New mass in solar units')
	    GSUR_NEW=GSUR_OLD*MASS_NEW/MASS_OLD
	    EXIT
	  END IF
	END DO
!
        VTURB=0.0D0; CALL GEN_IN(VTURB,'Turbulent velcity in km/s)')
        LOW_LIM=1; CALL GEN_IN(LOW_LIM,'Depth to begin revised mass estimate')
        HIGH_LIM=ND; CALL GEN_IN(HIGH_LIM,'Depth to end revised mass estimate')
!
        SUM_ERROR=0
        SUM_R=0
	DO I=1,ND
	  READ(STRING(I+1),*)R,V,E,VdVdR,dPdR,g_TOT,g_RAD,g_ELEC,Gamma
	  g_TOT=g_RAD-GSUR_NEW*(RND/R)**2
	  Gamma=g_RAD/GSUR_NEW*(R/RND)**2
	  IF(VTURB .EQ. 0)THEN
	    dTPdR=0.0D0
	  ELSE
	    dTPdR=-0.5D0*VTURB*VTURB*(2.0D0/R+VdVdR/V/V)
	  END IF
          DENOM=(ABS(VdVdR)+ ABS(dPdR)+ ABS(dTPdR)+ABS(g_TOT))
          E=200.0D0*(VdVdR+dPdR+dTPdR-g_TOT)/DENOM
          MT=g_rad/g_elec-1.0D0
          IF(I .GE. LOW_LIM .AND. I .LE. HIGH_LIM)THEN
            SUM_ERROR=SUM_ERROR+0.005*E/R**2/DENOM
            SUM_R=SUM_R+GRAV_CON/R**4/DENOM**2
          END IF
!
	  P_VEL(I)=V
	  P_REQ(I)=VdVdR+dPdR+dTPdR+GSUR_NEW*(RND/R)**2
	  P_GRAD(I)=g_RAD
	  P_GELEC(I)=g_ELEC
!
	  IF(R .GT. 9.99E+04)THEN
	    FMT='(1X,ES12.6,ES13.4,F9.2,6(ES14.4),2F11.2)'
	  ELSE
	    FMT='(1X,F12.6,ES13.4,F9.2,6(ES14.4),2F11.3)'
	  END IF                    
	  WRITE(LU_OUT,FMT)R,V,E,VdVdR,dPdR,dTPdR,g_TOT,g_RAD,g_ELEC,Gamma,MT
	END DO
!
        DEL_M=SUM_ERROR/SUM_R
        WRITE(T_OUT,*)'Better fit to hydrostatic equation will be obtained with a mass ',MASS_NEW-DEL_M
!
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Momentum equation is:',
	1                    ' VdV/dr = - dPdR/ROH - g + g_RAD'
	WRITE(LU_OUT,'(1X,A,A)')'        or          :',
	1                    ' VdV/dr = - dPdR/ROH + g_tot'
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Error is 200.0D0*(VdVdR+dPdR_ON_ROH-g_TOT)/',
	1        '( ABS(VdVdR)+ ABS(dPdR_ON_ROH)+ ABS(g_TOT) )'
C
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A)')'Gamma = g_rad/g [g=g_GRAV] '
C
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,ES14.4,A,F6.3,A)')
	1          'Surface gravity is: ',GSUR_NEW,' (i.e., log g=',LOG10(GSUR_NEW),')'
	WRITE(LU_OUT,'(1X,A,F8.2,A)')
	1          'Stars mass is: ',MASS_NEW,' Msun'
	CLOSE(LU_OUT)
!
	P_REQ(1:ND)=P_REQ(1:ND)/P_GELEC(1:ND)
	P_GRAD(1:ND)=P_GRAD(1:ND)/P_GELEC(1:ND)
!
	CALL DP_CURVE(ND,P_VEL,P_REQ)
	CALL DP_CURVE(ND,P_VEL,P_GRAD)
	CALL GRAMON_PGPLOT('V(km/s)','g/g\delec\u','(vdv/dr + \gr\u-1\d dP\dg\u/dr + g )/g\delec\u \p2',' ')
!
	DO I=1,ND
	  WRITE(40,'(I5,3ES16.6)')I,P_REQ(I),P_GELEC(I),P_GRAD(I)
	END DO
	CLOSE(UNIT=40)
	P_REQ(1:ND)=DLOG10(P_REQ(1:ND)*P_GELEC(1:ND))
	P_GRAD(1:ND)=DLOG10(P_GRAD(1:ND)*P_GELEC(1:ND))
	P_GELEC(1:ND)=DLOG10(P_GELEC(1:ND))
!
	CALL DP_CURVE(ND,P_VEL,P_REQ)
	CALL DP_CURVE(ND,P_VEL,P_GRAD)
	CALL DP_CURVE(ND,P_VEL,P_GELEC)
	CALL GRAMON_PGPLOT('V(km/s)','Log g','( vdv/dr + \gr\u-1\d dP/dr + g )',' ')
	STOP
	END
