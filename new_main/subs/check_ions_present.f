!
! Subroutine designed to facilitae (or check for inconsistencies) in
! the use of CMFGEN.
!
! Routine indicataes to the user whether ions can be deleted from the
! model, or whether it may be necesary to add additional high
! ionization stages.
!
	SUBROUTINE CHECK_IONS_PRESENT(ND)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created: 28-Jan-2009.
!
	INTEGER ND
	REAL*8 MAX_RATIO(NUM_IONS)
!
	REAL*8, PARAMETER ::  LOW_LIMIT=1.0D-10
	REAL*8, PARAMETER :: HIGH_LIMIT=1.0D-04
	REAL*8 T1,T2
	REAL*8 HDKT_EV
	REAL*8 T_VAL
	REAL*8 ED_VAL
	REAL*8 ION_FRAC
	REAL*8 INTERNAL_ENERGY
	REAL*8 EXTRA_INT_EN
	REAL*8 ATOM_ION_DENSITY
	INTEGER ISPEC
	INTEGER J,K
	INTEGER ID
	INTEGER DPTH_INDX
	LOGICAL FIRST_TIME
!
	INTEGER LUER,ERROR_LU,LUWARN,WARNING_LU
	EXTERNAL ERROR_LU,WARNING_LU
!
! The following table gives IP in eV. Grabbed off the web.
! As its only used to guide the user, the IP's need not be very accurate.
!
	REAL*8 IP(20,28)
	DATA IP(1:1,1)/    13.60/
	DATA IP(1:2,2)/    24.59,  54.42/
	DATA IP(1:3,3)/     5.39,  75.64, 122.45/
	DATA IP(1:4,4)/     9.32,  18.21, 153.90, 217.72/
	DATA IP(1:5,5)/     8.30,  25.16,  37.93, 259.37, 340.22/
	DATA IP(1:6,6)/    11.26,  24.38,  47.89,  64.49, 392.09, 489.99/
	DATA IP(1:7,7)/    14.53,  29.60,  47.45,  77.47,  97.89, 552.07, 667.04/
	DATA IP(1:8,8)/    13.62,  35.12,  54.94,  77.41, 113.90, 138.12, 739.28, 871.41/
	DATA IP(1:9,9)/    17.42,  34.97,  62.71,  87.14, 114.24, 157.16, 185.19, 953.91,1103.11/
	DATA IP(1:10,10)/  21.56,  40.96,  63.45,  97.12, 126.21, 157.93, 207.27, 239.10,1195.82,1362.20/
	DATA IP(1:11,11)/   5.14,  47.28,  71.62,  98.91, 138.40, 172.18, 208.50, 264.25, 299.86,1465.11,1648.71/
	DATA IP(1:12,12)/   7.65,  15.04,  80.14, 109.27, 141.26, 186.76, 225.02, 265.96, 328.06, 367.50,1761.80,
	1                1962.66/
	DATA IP(1:13,13)/   5.99,  18.83,  28.45, 119.99, 153.83, 190.48, 241.76, 284.65, 330.13, 398.74, 442.00,
	1                2315.02,2304.14/
	DATA IP(1:14,14)/   8.15,  16.35,  33.49,  45.14, 166.77, 205.26, 246.46, 303.54, 351.12, 401.37, 476.36,
	1                 523.42,2437.63,2673.18/
	DATA IP(1:15,15)/  10.49,  19.76,  30.20,  51.44,  65.02, 220.42, 263.57, 309.60, 372.13, 424.42, 479.46,
	1                 560.81, 611.74,2816.91,3069.84/
	DATA IP(1:16,16)/  10.36,  23.34,  34.79,  47.22,  72.59,  88.05, 280.94, 328.74, 379.55, 447.50, 504.84,
	1                 564.44, 652.22, 707.01,3223.78,3494.19/
	DATA IP(1:17,17)/  12.97,  23.82,  39.61,  53.47,  67.80,  97.03, 114.19, 348.28, 400.06, 455.62, 529.28,
	1                 592.00, 656.71, 749.76, 809.40,3658.52,3946.30/
	DATA IP(1:18,18)/  15.76,  27.63,  40.74,  59.81,  75.02,  91.01, 124.32, 143.46, 422.45, 478.68, 538.96,
	1                 618.26, 686.10, 755.74, 854.77, 918.03,4120.88,4426.22/
	DATA IP(1:19,19)/   4.34,  31.63,  45.81,  60.91,  82.66,  99.39, 117.56, 154.88, 175.82, 503.81, 564.75,
	1                 629.42, 714.62, 786.65, 861.06, 968.02,1033.42,4610.85,4934.04/
	DATA IP(1:20,20)/   6.11,  11.87,  50.91,  67.27,  84.50, 108.78, 127.17, 147.23, 188.54, 211.28, 591.90,
	1                 657.20, 726.64, 817.64, 894.54, 974.24,1087.21,1157.80,5128.76,5469.86/
	DATA IP(1:20,21)/   6.56,  12.80,  24.76,  73.49,  91.65, 110.68, 137.95, 158.06, 180.03, 225.17, 249.80,
	1                 687.36, 756.69, 830.80, 927.50,1009.48,1094.47,1212.62,1287.97,5674.75/
	DATA IP(1:20,22)/   6.83,  13.58,  27.49,  43.27,  99.30, 119.53, 140.85, 170.39, 192.05, 215.92, 265.07,
	1                 291.49, 787.84, 863.14, 941.90,1043.68,1130.74,1220.91,1346.32,1425.40/
	DATA IP(1:20,23)/   6.75,  14.66,  29.33,  46.71,  65.28, 128.13, 150.59, 173.39, 205.83, 230.50, 255.69,
	1                 308.13, 336.28, 895.99, 976.00,1060.26,1168.05,1260.29,1354.61,1486.24/
	DATA IP(1:20,24)/   6.77,  16.49,  30.96,  49.16,  69.46,  90.63, 160.18, 184.69, 209.25, 244.39, 270.82,
	1                 297.97, 354.77, 384.16,1010.62,1096.54,1184.64,1298.64,1396.07,1495.56/
	DATA IP(1:20,25)/   7.43,  15.64,  33.66,  51.20,  72.45,  95.56, 119.19, 194.54, 221.80, 248.33, 285.95,
	1                 314.35, 343.58, 402.96, 435.16,1134.68,1224.02,1317.30,1436.49,1539.09/
	DATA IP(1:20,26)/   7.90,  16.19,  30.65,  54.83,  75.04,  99.08, 124.99, 151.11, 233.61, 262.11, 290.20,
	1                 330.83, 360.99, 392.18, 457.06, 489.26,1266.51,1357.72,1456.18,1581.59/
	DATA IP(1:20,27)/   7.88,  17.08,  33.50,  51.30,  79.49, 101.98, 128.93, 157.85, 186.13, 275.38, 304.71,
	1                 335.80, 379.33, 411.46, 443.59, 511.95, 546.58,1397.21,1504.58,1603.35/
	DATA IP(1:20,28)/   7.64,  18.17,  35.19,  54.93,  76.06, 107.79, 132.66, 161.68, 192.78, 224.59, 320.98,
	1                 352.38, 384.51, 430.12, 464.32, 498.52, 571.08, 607.03,1541.17,1647.92/
!
	MAX_RATIO=1.0D0
	LUER=ERROR_LU()
	LUWARN=WARNING_LU( )
!
! Determine ionzation fractions.
!
	DO ISPEC=1,NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      T1=0.0D0
	      DO K=1,ND
	        T2=SUM(ATM(ID)%XzV_F(:,K))
	        T1=MAX(T1,T2/POP_SPECIES(K,ISPEC))
	      END DO
	      MAX_RATIO(ID)=T1
	    END DO
	    ID=SPECIES_END_ID(ISPEC)-1
	    MAX_RATIO(ID+1)=MAXVAL(ATM(ID)%DXzV_F/POP_SPECIES(:,ISPEC))	   
	  END IF
	END DO
!
! Check whether the lowest ionization stages might be omitted.
!
	WRITE(LUWARN,'(/,/,A,A)')' Checking whether some low ionization stages ',
	1                        'may be omitted from the model.'
	DO ISPEC=1,NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	      IF(MAX_RATIO(ID) .GT. LOW_LIMIT)THEN
	        DO K=SPECIES_BEG_ID(ISPEC), ID-1
	          WRITE(LUWARN,'(3X,A,T11,A,ES8.2,A,ES8.2)')TRIM(ION_ID(K)),
	1           ' may not need to be include in model. Its maximum fractional abundance of ',
	1           MAX_RATIO(K),' < ',LOW_LIMIT
	        END DO
	        EXIT
	      END IF
	    END DO
	  END IF
	END DO
!
! Check whether the highest ionization stages might be omitted.
!
	WRITE(LUWARN,'(A,A)')' Checking whether some high ionization stages ',
	1                        'may be omitted from the model.'
	DO ISPEC=1,NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	      IF(MAX_RATIO(ID) .GT. HIGH_LIMIT)THEN
	        DO K=ID+1,SPECIES_END_ID(ISPEC)
	          IF(K .EQ. SPECIES_END_ID(ISPEC))THEN
	            J=NINT(ATM(K-1)%ZXzV)+1
	          ELSE
	            J=NINT(ATM(K)%ZXzV)
	          END IF
	          WRITE(LUWARN,'(3X,A,T11,A,ES8.2,A,ES8.2)')TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J)),
	1             'may not need to be include in model. Its maximum fractional abundance is ',
	1            MAX_RATIO(K),' < ',HIGH_LIMIT
	        END DO
	        EXIT
	      END IF
	    END DO
	  END IF
	END DO
!
! Check to see if lower ionization species need to be included.
!
	FIRST_TIME=.TRUE.
	DO ISPEC=1,NUM_SPECIES
	  ID=SPECIES_BEG_ID(ISPEC)
	  IF(SPECIES_PRES(ISPEC) .AND. ATM(ID)%ZXzV .GT. 1.1D0)THEN
	    IF(FIRST_TIME)THEN
	      WRITE(LUWARN,'(A,A)')' Checking whether lower ionization stages ',
	1                             'need to be included in model'
	      FIRST_TIME=.FALSE.
	    END IF
	    IF(MAX_RATIO(ID) .GT. 0.1D0)THEN
	         K=NINT(ATM(ID)%ZXzV)-1
	         WRITE(LUWARN,'(3X,A,T11,A,A,A,ES8.1,A)')TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(K)),
	1              'may need to be included in the model (Maximum ',TRIM(ION_ID(ID)),
	1              ' ionization fraction is',MAX_RATIO(ID),')'
	    ELSE
	         EXIT
	    END IF
	  END IF
	END DO
!
! Now do a quick and dirty check to see whether additional ionization stages should be added.
! We first get the maximum temperature (in SN models this may not be at the inner
! boundary).
!
	T_VAL=0.0D0
	DO K=1,ND
	  IF(T(K) .GT. T_VAL)THEN
	    T_VAL=T(K)
	    ED_VAL=ED(K)
	    DPTH_INDX=K
	  END IF
	END DO
!
! To check whether additional ionization stages need to be included,
! we allow for the possibility that the previous ionization stage
! might not be dominant. Thus we multipy by ION_FRAC when its less
! than 1.
! 
	HDKT_EV=4.7994145D0/4.13566733D0
	FIRST_TIME=.TRUE.
	EXTRA_INT_EN=0.0D0
	DO ISPEC=1,NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    ID=SPECIES_END_ID(ISPEC)-1
	    J=ATM(ID)%ZXzV+2
	    K=NINT(AT_NO(ISPEC))
	    ION_FRAC=1.0D0
	    DO ID=J,MIN(K,20)
	       ION_FRAC=ION_FRAC*(T_VAL**1.5D0)*EXP(-HDKT_EV*IP(ID,K)/T_VAL)/2.07D-22/ED_VAL
	       IF(ION_FRAC .GT. HIGH_LIMIT)THEN
	         IF(FIRST_TIME)THEN
	           WRITE(LUWARN,'(A)')' '
	           WRITE(LUWARN,'(A,A)')' Checking whether additional higher ionization stages ',
	1                             ' need to be included in the model.'
	           WRITE(LUWARN,'(A,A,/,A)')' NB: If XzV needs to be included it will also be necessary to',
	1                             ' include XzIV, as this was',
	1                             '      only included as the ground state.'
	           WRITE(LUWARN,'(A,1X,A,I3,A,ES10.4,A,ES10.4)')' Parameters at check depth:',
	1	          'Depth=',DPTH_INDX,'   T=',T_VAL,'   ED=',ED_VAL
	           FIRST_TIME=.FALSE.
	         END IF
	         WRITE(LUWARN,'(3X,A,T11,A,ES8.1,A)')TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(ID)),
	1              'may need to be included in the model (IF ~',ION_FRAC,')'
	         IF(ION_FRAC .GT. 1.0D0)ION_FRAC=1.0D0
	         EXTRA_INT_EN=EXTRA_INT_EN+ION_FRAC*POP_SPECIES(DPTH_INDX,ISPEC)*IP(ID,K)
	       ELSE
	         EXIT
	       END IF
	    END DO
	  END IF
	END DO
!
	ATOM_ION_DENSITY=0.0D0
	DO ISPEC=1,NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    ATOM_ION_DENSITY=ATOM_ION_DENSITY+POP_SPECIES(DPTH_INDX,ISPEC)
	  END IF
	END DO
!
	INTERNAL_ENERGY=0.0D0
	DO ISPEC=1,NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      T1=0.0D0
	      K=DPTH_INDX
	      T2=SUM(ATM(ID)%XzV_F(:,K))
	      J=ATM(ID)%ZXzV-1
	      IF(J .NE. 0)INTERNAL_ENERGY=INTERNAL_ENERGY+T2*SUM(IP(1:J,NINT(AT_NO(ISPEC))))
!	      WRITE(6,*)ISPEC,T2/ATOM_ION_DENSITY,SUM(IP(1:J,NINT(AT_NO(ISPEC))))
	    END DO
	    J=J+1
	    ID=SPECIES_END_ID(ISPEC)-1
!	    WRITE(6,*)ISPEC,ATM(ID)%DXzV_F(DPTH_INDX)/ATOM_ION_DENSITY,SUM(IP(1:J,NINT(AT_NO(ISPEC))))
	    INTERNAL_ENERGY=INTERNAL_ENERGY+ATM(ID)%DXzV_F(DPTH_INDX)*SUM(IP(1:J,NINT(AT_NO(ISPEC))))
	  END IF
	END DO
!
	T1=1.5D0*0.86173D0*T(DPTH_INDX)*(ED(DPTH_INDX)+ATOM_ION_DENSITY)
	WRITE(LUWARN,'(A)')' '
	WRITE(LUWARN,'(A,ES9.3,A)')'                   Temperature (10^4 K) is ',T(DPTH_INDX)
	WRITE(LUWARN,'(A,ES9.3,A)')'           Atom/ion density (per cm^3)) is ',ATOM_ION_DENSITY
	WRITE(LUWARN,'(A,ES9.3,A)')'                 Electrons per atom/ion is ',ED(DPTH_INDX)/ATOM_ION_DENSITY
	WRITE(LUWARN,'(A,ES9.3,A)')' Total thermal kinetic energy (ev/atom) is ',T1/ATOM_ION_DENSITY
	WRITE(LUWARN,'(A,ES9.3,A)')'  Approximate internal energy (ev/atom) is ',INTERNAL_ENERGY/ATOM_ION_DENSITY
	WRITE(LUWARN,'(A,ES9.3,A)')' App. missing internal energy (ev/atom) is ',EXTRA_INT_EN/ATOM_ION_DENSITY
!
! 4/c . sigma T^4 conerted to eV  (NB: Int J = sigma T^4 / pi).
!
	T1=4.0D0*3.14159265D0*5.670400D-05*1.0D+16/1.60217733D-12/2.99792458D+10
	WRITE(LUWARN,'(A,ES9.3,A)')'             Radiation energy (ev/atom) is ',(T1*T(DPTH_INDX)**4)/ATOM_ION_DENSITY
	FLUSH(UNIT=6)
!
	RETURN
	END
