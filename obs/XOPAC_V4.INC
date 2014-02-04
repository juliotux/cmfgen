C
C NB: This differs from XOPAC.RAW used by CMFGEN in that the full
C arrays, not the super level atoms, are accessed.
C
C Compute opacities and emissivities due to K shell ionization. In all cases
C it is assumed that 2 electrons are ejected.
C
C Altered: Bux fix. NCI replaced by NXzV. Provided CI was included in
c                     Compile code, minimal effect since opacity dominated by
C                     lowest levels.
C
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	      T2=AT_NO(SPECIES_LNK(ID))+1-ATM(ID)%ZXzV
!	      T1=XCROSS(Fl,AT_NO(SPECIES_LNK(ID)),T2)
	      T1=XCROSS_V2(Fl,AT_NO(SPECIES_LNK(ID)),T2,IZERO,IZERO,L_FALSE,L_FALSE)
	      IF(T1 .NE. 0.0D0)THEN
	        DO I=1,ND
	          T2=0.0D0			!Temporary CHI
	          T3=0.0D0			!Temporary ETA
	          T4=ATM(ID+1)%XzVLTE_F(1,I)/ATM(ID+1)%XzV_F(1,I)*EMHNUKT(I)
	          DO J=1,ATM(ID)%NXzV_F
		    T2=T2+(ATM(ID)%XzV_F(J,I)-ATM(ID)%XzVLTE_F(J,I)*T4)
	            T3=T3+ATM(ID)%XzVLTE_F(J,I)
	          END DO
	          CHI(I)=CHI(I)+T1*T2
	          ETA(I)=ETA(I)+T1*T3*T4*TWOHCSQ*(FL**3)
	        END DO
	     END IF
	    END IF
	END DO
C
