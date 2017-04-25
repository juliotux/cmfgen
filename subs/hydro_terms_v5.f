!
! Subroutine to compare terms in the Hydronamical equations. Routine is
! meant for illustration purposes only.
!
	SUBROUTINE HYDRO_TERMS_V5(POP_ATOM,R,V,T,SIGMA,ED,CLUMP_FAC,LUM_STAR,
	1              LOGG,STARS_MASS,MEAN_ATOMIC_WEIGHT,
	1              FLUX_MEAN_OPAC,ROSS_MEAN_OPAC,ELEC_MEAN_OPAC,
	1              LAM_FLUXMEAN_BAND_END,BAND_FLUXMEAN,
	1              PRESSURE_VTURB,PLANE_PARALLEL,PLANE_PARALLEL_NOV,
	1              LST_ITERATION,BAND_FLUX,N_FLUXMEAN_BANDS,LU_OUT,ND)
	IMPLICIT NONE
!
! Altered 01-Jan-2015 : Turbulent pressure term added, and also output when non-zero.
! Altered 26-Jan-2014 : Changed ouput to TYPE_ATM and increased its length.
! Altered 06-Jan-2014 : Changed to V5 (inserted LOGG, ROSS_MEAN_OPAC and CLUMP_FAC in call).
! Altered 31-Mar-2010 : Changed from V3 to V4: Inserted LST_ITERATION in call.
! Altered 11-Mar-2008 : Changed from V2 to V3
!                       PRESSURE_VTURB,PLANE_PARALLEL & PLANE_PARALLEL_NOV inserted into call.
! Altered 16-Jan-2005 : Accuracy of dPdR/ROH improved. Now use monotonic
!                         interpolation. Depth index added output.
! Altered 17-Mar-2001 : Altered V output format
! Altered 27-Oct-2000 : Percentage error improved.
!                       STARS mass now output.
! Altered 24-May-1996 : L GEN_ASCI_OPEN and ERROR_LU installed.
! Altered 06-Mar-1996 : Description of terms added to output.
! Created 28-Mar-1994
!
	INTEGER ND,LU_OUT
!
	REAL*8 POP_ATOM(ND)
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 T(ND)
	REAL*8 SIGMA(ND)
	REAL*8 ED(ND)
	REAL*8 CLUMP_FAC(ND)
	REAL*8 LUM_STAR(ND)                        
!
	REAL*8 FLUX_MEAN_OPAC(ND)
	REAL*8 ROSS_MEAN_OPAC(ND)
	REAL*8 ELEC_MEAN_OPAC(ND)
!
	INTEGER N_FLUXMEAN_BANDS
	REAL*8 LAM_FLUXMEAN_BAND_END(ND)
	REAL*8 BAND_FLUXMEAN(ND,N_FLUXMEAN_BANDS)
	REAL*8 BAND_FLUX(ND,N_FLUXMEAN_BANDS)
!
	REAL*8 LOGG
	REAL*8 STARS_MASS
	REAL*8 MEAN_ATOMIC_WEIGHT
	REAL*8 PRESSURE_VTURB
	LOGICAL PLANE_PARALLEL
	LOGICAL PLANE_PARALLEL_NOV
	LOGICAL LST_ITERATION 
!
! External functions
!
	REAL*8 BOLTZMANN_CONSTANT
	REAL*8 GRAVITATIONAL_CONSTANT
	REAL*8 MASS_SUN
	REAL*8 LUM_SUN
	REAL*8 FUN_PI
	REAL*8 SPEED_OF_LIGHT
	REAL*8 ATOMIC_MASS_UNIT
!
	EXTERNAL BOLTZMANN_CONSTANT,GRAVITATIONAL_CONSTANT,MASS_SUN
	EXTERNAL LUM_SUN,FUN_PI,SPEED_OF_LIGHT,ATOMIC_MASS_UNIT
!
	REAL*8 MOD_PRESSURE(ND)
	REAL*8 COEF(ND,4)
	REAL*8 TAU(ND)
	REAL*8 OPAC(ND)
	REAL*8 TB(ND)
	REAL*8 TC(ND)
!
! Local variables.
!
	REAL*8 VdVdR
	REAL*8 dPdR_ON_ROH
	REAL*8 dPTURBdR_ON_ROH
	REAL*8 g_TOT
	REAL*8 g_ELEC
	REAL*8 g_GRAV
	REAL*8 g_RAD
	REAL*8 ERROR
	REAL*8 RLUMST(ND)
!
	REAL*8 GRAV_CON
	REAL*8 dP_CON
	REAL*8 PTURB_CON
	REAL*8 RAD_CON
	REAL*8 T1
	REAL*8 RPHOT
	REAL*8 GPHOT
	REAL*8 ERROR_MAX
	REAL*8 ERROR_SUM
	REAL*8 ERROR_SQ
!
	CHARACTER(LEN=12) TYPE_ATM
	CHARACTER(LEN=80) FMT
	INTEGER I,J,IOS,ERROR_LU
	INTEGER ERROR_CNT
	EXTERNAL ERROR_LU
	INTEGER, PARAMETER :: IONE=1
!
	I=MAX(132,(N_FLUXMEAN_BANDS+1)*12+16)
	CALL GEN_ASCI_OPEN(LU_OUT,'HYDRO','UNKNOWN',' ',' ',I,IOS)
	IF(IOS .NE. 0)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error opening HYDRO in HYDRO_TERMS'
	  WRITE(I,*)'IOS=',IOS
	END IF
!
	ERROR_SUM=0.0D0
	ERROR_SQ=0.0D0
	ERROR_MAX=0.0D0
	ERROR_CNT=0
!
! 1.0D-06 = [ 10^(-10) {from 1/r) * 10^4 {from T}
!
	dP_CON=1.0D-06*BOLTZMANN_CONSTANT()/MEAN_ATOMIC_WEIGHT/
	1               ATOMIC_MASS_UNIT()
	PTURB_CON=1.0D+10*0.5D0*PRESSURE_VTURB*PRESSURE_VTURB
!
	GRAV_CON=1.0D-20*GRAVITATIONAL_CONSTANT()*STARS_MASS*MASS_SUN()
!                               
!       1.0E-10 (CHI) / 1.0E-20 (1/r^2)
!
	RAD_CON=1.0D-30*LUM_SUN()/ATOMIC_MASS_UNIT()/4.0D0/FUN_PI()/
	1                  SPEED_OF_LIGHT()/MEAN_ATOMIC_WEIGHT
!
	TYPE_ATM='P'
        T1=LOG(POP_ATOM(5)/POP_ATOM(1))/LOG(R(1)/R(5))
        WRITE(TYPE_ATM(2:),'(ES10.3)')T1
!
	DO I=1,ND
	  OPAC(I)=MAX(ROSS_MEAN_OPAC(I),ELEC_MEAN_OPAC(I))*CLUMP_FAC(I)
	END DO
	CALL TORSCL(TAU,OPAC,R,TB,TC,ND,'LOGMON',TYPE_ATM)
	T1=0.6667D0
	CALL MON_INTERP(RPHOT,IONE,IONE,T1,IONE,R,ND,TAU,ND)
	IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NOV)RPHOT=R(ND)
	GPHOT=GRAV_CON/RPHOT/RPHOT	
!                                              
! Output file header.
!
	IF(PRESSURE_VTURB .EQ. 0.0D0)THEN
	  WRITE(LU_OUT,'(1X,6X,A,6X, 8X,A,3X, 2X,A, 5(4X,A,1X), 4X,A,2X,A)')
	1     'R','V','% Error','    VdVdR',
	1                       ' dPdR/ROH',
	1                       '    g_TOT',
	1                       '    g_RAD',
	1                       '   g_ELEC','Gamma','Depth'
	ELSE
	  WRITE(LU_OUT,'(1X,6X,A,6X, 8X,A,3X, 2X,A, 6(4X,A,1X), 4X,A,2X,A)')
	1     'R','V','% Error','    VdVdR',
	1                       ' dPdR/ROH',
	1                       'dTPdR/ROH',
	1                       '    g_TOT',
	1                       '    g_RAD',
	1                       '   g_ELEC','Gamma','Depth'
	END IF
! 
	DO I=1,ND
	  MOD_PRESSURE(I)=(POP_ATOM(I)+ED(I))*T(I)
	END DO
	CALL MON_INT_FUNS_V2(COEF,MOD_PRESSURE,R,ND)
!	                                     
	DO I=1,ND
	  VdVdR=V(I)*V(I)*(SIGMA(I)+1.0D0)/R(I)
	  dPdR_ON_ROH=dP_CON*COEF(I,3)/POP_ATOM(I)
	  dPTURBdR_ON_ROH=PTURB_CON*1.0D-10*(-2.0D0/R(I)-(SIGMA(I)+1.0D0)/R(I))
!
! Gravitational force per unit mass (in cgs units) and radiation pressure forces.
!
	  IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NOV)THEN
	    g_grav=GRAV_CON/R(ND)/R(ND)
	    g_rad=RAD_CON*LUM_STAR(ND)*FLUX_MEAN_OPAC(I)/POP_ATOM(I)/R(ND)/R(ND)
	    g_elec=RAD_CON*LUM_STAR(ND)*ELEC_MEAN_OPAC(I)/POP_ATOM(I)/R(ND)/R(ND)
	  ELSE
	    g_grav=GRAV_CON/R(I)/R(I)
	    g_rad=RAD_CON*LUM_STAR(ND)*FLUX_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
	    g_elec=RAD_CON*LUM_STAR(ND)*ELEC_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
!	    g_rad=RAD_CON*LUM_STAR(I)*FLUX_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
!	    g_elec=RAD_CON*LUM_STAR(I)*ELEC_MEAN_OPAC(I)/POP_ATOM(I)/R(I)/R(I)
	  END IF
!
	  g_TOT= g_RAD-g_GRAV
	  ERROR=200.0D0*(VdVdR+dPdR_ON_ROH+dPTURBdR_ON_ROH-g_TOT)/
	1        ( ABS(VdVdR)+ ABS(dPdR_ON_ROH)+ ABS(dPTURBdR_ON_ROH)+ ABS(g_TOT) )
!
	  IF(V(I) .LT. 5.0D0)THEN
	    ERROR_SUM=ERROR_SUM+ERROR
	    ERROR_SQ=ERROR_SQ+ERROR*ERROR
	    ERROR_MAX=MAX(ERROR_MAX,ABS(ERROR))
	    ERROR_CNT=ERROR_CNT+1
	  END IF
!
	  IF(PRESSURE_VTURB .EQ. 0.0D0)THEN
	    IF(R(I) .GT. 9.99D+04)THEN
	      FMT='(1X,ES12.6,ES13.4,F9.2,5(ES14.4),F9.2,I7)'
	    ELSE
	      FMT='(1X,F12.6,ES13.4,F9.2,5(ES14.4),F9.2,I7)'
	    END IF
	    WRITE(LU_OUT,FMT)
	1             R(I),V(I),ERROR,VdVdR,dPdR_ON_ROH,
	1             g_TOT,g_RAD,g_ELEC,g_RAD/g_GRAV,I
	  ELSE                  
	    IF(R(I) .GT. 9.99D+04)THEN
	      FMT='(1X,ES12.6,ES13.4,F9.2,6(ES14.4),F9.2,I7)'
	    ELSE
	      FMT='(1X,F12.6,ES13.4,F9.2,6(ES14.4),F9.2,I7)'
	    END IF
	    WRITE(LU_OUT,FMT)
	1             R(I),V(I),ERROR,VdVdR,dPdR_ON_ROH,dPTURBdR_ON_ROH,
	1             g_TOT,g_RAD,g_ELEC,g_RAD/g_GRAV,I
	  END IF                    
!
	END DO
!
! Write out a summary of momentum euqation etc to remind the user of the
! meaning of individual terms.
!
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Momentum equation is:',
	1                    ' VdV/dr = - dPdR/ROH - dTPdR/ROH - g + g_RAD'
	WRITE(LU_OUT,'(1X,A,A)')'        or          :',
	1                    ' VdV/dr = - dPdR/ROH - dTPdR/ROH + g_tot'
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Error is 200.0D0*(VdVdR+dPdR_ON_ROH+dTPdR/ROH-g_TOT)/',
	1        '( ABS(VdVdR)+ ABS(dPdR_ON_ROH)+ ABS(dTPdR/ROH)+ ABS(g_TOT) )'
!
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A)')'Gamma = g_rad/g [g=g_GRAV] '
!
	WRITE(LU_OUT,'(1X,A)')' '
	T1=GRAV_CON/RPHOT/RPHOT
	WRITE(LU_OUT,'(1X,A,ES14.4,A,4X,ES10.4,A)')
	1          '         Photospheric radius is: ',RPHOT,'(10^10 cm)',RPHOT/6.96D0,'(Rsun)'
	WRITE(LU_OUT,'(1X,A,ES14.4,A,F7.4,A)')
	1          'Photospheric surface gravity is: ',T1,' (',LOG10(T1),')'
	WRITE(LU_OUT,'(1X,A,ES14.4,A,F7.4,A)')
	1          '   Specified surface gravity is: ',10**LOGG,' (',LOGG,')'
	WRITE(LU_OUT,'(1X,A,F8.2,A)')
	1          '                  Stars mass is: ',STARS_MASS,' Msun'
	WRITE(LU_OUT,'(A,A)')' Power law exponent at outer boundary (TYPE_ATM)=',TRIM(TYPE_ATM(2:))

!
	RLUMST(:)=BAND_FLUX(:,N_FLUXMEAN_BANDS)
        DO J=N_FLUXMEAN_BANDS,2,-1
          BAND_FLUXMEAN(:,J)=BAND_FLUXMEAN(:,J)-BAND_FLUXMEAN(:,J-1)
          BAND_FLUX(:,J)=BAND_FLUX(:,J)-BAND_FLUX(:,J-1)
        END DO
        DO J=1,N_FLUXMEAN_BANDS
          BAND_FLUXMEAN(:,J)=BAND_FLUXMEAN(:,J)/FLUX_MEAN_OPAC(:)
          BAND_FLUX(:,J)=BAND_FLUX(:,J)/RLUMST(:)
        END DO
!
	WRITE(LU_OUT,'(A)')CHAR(12)                         !Formfeed
        WRITE(LU_OUT,'(3X,A,5X,A,6X,14F12.2)')'I','V',
	1                      (LAM_FLUXMEAN_BAND_END(J),J=1,N_FLUXMEAN_BANDS)
        DO I=1,ND
          WRITE(LU_OUT,'(X,I3,15ES12.3)')I,V(I),(BAND_FLUXMEAN(I,J),J=1,N_FLUXMEAN_BANDS)
          WRITE(LU_OUT,'(16X,15ES12.3)')(BAND_FLUX(I,J),J=1,N_FLUXMEAN_BANDS)
	  T1=FLUX_MEAN_OPAC(I)/ELEC_MEAN_OPAC(I)
	  DO J=1,N_FLUXMEAN_BANDS
	    IF(ABS(BAND_FLUX(I,J)) .LT. 1.0D-80)BAND_FLUX(I,J)=1.0D-80
          END DO 
	  WRITE(LU_OUT,'(16X,15ES12.3)')(T1*BAND_FLUXMEAN(I,J)/BAND_FLUX(I,J),J=1,N_FLUXMEAN_BANDS)
        END DO
!
	CLOSE(LU_OUT)
!
	IF(LST_ITERATION)THEN
	  ERROR_SQ=SQRT(ERROR_SQ/MAX(1,ERROR_CNT))
	  ERROR_SUM=ERROR_SUM/MAX(1,ERROR_CNT)
	  T1=ABS(LOGG-LOG10(GPHOT))
	  IF(ERROR_SQ .GT. 5.0D0 .OR. ERROR_MAX .GT. 20.0D0 .OR. T1 .GT. 0.02D0)THEN
	    I=ERROR_LU()
	    WRITE(I,*)' '
	    WRITE(I,*)'*****************************************************************************************'
	    WRITE(I,'(A)')' Possible error with hydrostatic structure --- large error in photosphere'
	    WRITE(I,'(A,ES10.2)')'              Mean error is',ERROR_SUM
	    WRITE(I,'(A,ES10.2)')' Root mean squared error is',ERROR_SQ
	    WRITE(I,'(A,ES10.2)')'           Maximum error is',ERROR_MAX
	    WRITE(I,'(A,ES10.2)')'       Stars mass (Msun) is',STARS_MASS
	    WRITE(I,'(A,ES10.2)')'        R(phot, 10^10cm) is',RPHOT
	    WRITE(I,'(A,ES10.2)')'             log G(phot) is',LOG10(GPHOT)
	    WRITE(I,'(A,ES10.2)')'         Specified log g is',LOGG
	    WRITE(I,*)'*****************************************************************************************'
	    WRITE(I,*)' '
	  END IF
	END IF
!
	RETURN
	END
