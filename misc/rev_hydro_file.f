!
! Auxilary program designed to modify the HYDRO file output from CMFGEN.
! Progam modifies the adopted stellar mass. Not that the percentage error
! os now defined so that it has a rang of pm 200%.
!
	PROGRAM REV_HYDRO_FILE
	USE GEN_IN_INTERFACE
!
! Cleaned: 07-Nov-200
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
	REAL*8 GRAVITATIONAL_CONSTANT, MASS_SUN
	EXTERNAL GRAVITATIONAL_CONSTANT, MASS_SUN
!
	CHARACTER*132 STRING(500)
	CHARACTER*132 FMT
	CHARACTER*132 FILENAME
!
	GRAV_CON=1.0D-20*GRAVITATIONAL_CONSTANT()*MASS_SUN()
!
	FILENAME='HYDRO'
	CALL GEN_IN(FILENAME,'Input hydro file')
	OPEN(UNIT=10,FILE='HYDRO',ACTION='READ',STATUS='OLD')
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
        I=0
!
	FILENAME='REV_HYDRO'
	CALL GEN_IN(FILENAME,'Output hydro file')
	CALL GEN_ASCI_OPEN(LU_OUT,FILENAME,'UNKNOWN',' ',' ',I,IOS)
	WRITE(LU_OUT,'(1X,4X,A,5X, 8X,A,3X, 2X,A, 5(5X,A,1X), 6X,A, 6X,A)')
	1        'R','V','% Error','   VdVdR',
	1                   'dPdR/ROH',
	1                   '   g_TOT',
	1                   '   g_RAD',
	1                    '  g_ELEC','Gamma','M(t)'
!
	READ(STRING(ND+1),*)RND
	DO I=1,NSTR
	  IF(INDEX(STRING(I),'Surface gravity is:') .NE. 0)THEN
	    J=INDEX(STRING(I),':')
	    READ(STRING(I)(J+1:),*)GSUR_OLD
	    MASS_OLD=GSUR_OLD/GRAV_CON*RND*RND
	    WRITE(T_OUT,'(1X,A,F6.2)')'Old mass is ',MASS_OLD
	    MASS_NEW=MASS_OLD
	    CALL GEN_IN(MASS_NEW,'New mass in solar units')
	    GSUR_NEW=GSUR_OLD*MASS_NEW/MASS_OLD
	    EXIT
	  END IF
	END DO
!
	DO I=1,ND
	  READ(STRING(I+1),*)R,V,E,VdVdR,dPdR,g_TOT,g_RAD,g_ELEC,Gamma
	  Gamma=Gamma*MASS_OLD/MASS_NEW
	  g_tot=(Gamma-1.0D0)*g_rad/Gamma
	  E=200.0D0*(VdVdR+dPdR-g_TOT)/(ABS(VdVdR)+ ABS(dPdR)+ ABS(g_TOT))
	  MT=g_rad/g_elec-1.0D0
!
	  IF(R .GT. 9.99E+04)THEN
	    FMT='(1X,1PE10.4,0PF12.2,F9.2,1P,5(E14.4),0P,2F11.2)'
	  ELSE
	    FMT='(1X,F10.4,F12.2,F9.2,1P,5(E14.4),0P,2F11.3)'
	  END IF                    
	  WRITE(LU_OUT,FMT)R,V,E,VdVdR,dPdR,g_TOT,g_RAD,g_ELEC,Gamma,MT
	END DO
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
	WRITE(LU_OUT,'(1X,A,1PE14.4)')
	1          'Surface gravity is: ',GSUR_NEW
	WRITE(LU_OUT,'(1X,A,F8.2,A)')
	1          'Stars mass is: ',MASS_NEW,' Msun'
!
	STOP
	END
