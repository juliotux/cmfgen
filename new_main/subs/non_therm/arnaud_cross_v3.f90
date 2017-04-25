	SUBROUTINE ARNAUD_CROSS_V3()
	USE MOD_NON_THERM
	USE MOD_CMFGEN
	IMPLICIT NONE
!
!	INTEGER NKT
!	REAL*8 XKT(NKT)
!
	INTEGER IT
	INTEGER IKT
	INTEGER IKT_ST
	INTEGER IKT_END
	REAL*8 U1
	REAL*8 T1,T2,T3
	INTEGER ID
	INTEGER I,J
!	INTEGER MAX_NUM_IONS
	LOGICAL NEGATIVE_CROSS_SEC
!	CHARACTER*12 ION_ID(MAX_NUM_IONS)
!
! Fitting formula from Arnaud & Rothenflug 1985:
! http://adsabs.harvard.edu/abs/1985A%26AS...60..425A
! Gives the cross section for direct ionisation.
!
! Note: The AR constants are in units of 1.0E-14 cm^2 eV^2.
!
	DO IT=1,NUM_THD
	  IF(THD(IT)%PRES)THEN
	    NEGATIVE_CROSS_SEC=.FALSE.
	    THD(IT)%CROSS_SEC=0.0D0
	    ID=THD(IT)%LNK_TO_ION
!
	    IF(THD(IT)%NTAB .EQ. 0)THEN
	      DO IKT=1,NKT
	        U1 = XKT(IKT) / THD(IT)%ION_POT
	        IF (U1 .GT. 1.0D0) THEN
	          T1=(1.0D0-1.0D0/U1)
	          T2=DLOG(U1)
	          T3 = 1.0D-14 * ( THD(IT)%A_COL*T1 + THD(IT)%B_COL*T1*T1 + &
	              THD(IT)%C_COL*T2 + THD(IT)%D_COL*T2/U1  ) / U1 / THD(IT)%ION_POT**2
	          IF(T3 .GE. 0.0D0)THEN
	            THD(IT)%CROSS_SEC(IKT)=T3
	          ELSE
	            THD(IT)%CROSS_SEC(IKT)=0.0D0
	            NEGATIVE_CROSS_SEC=.TRUE.
	          END IF
	        END IF
	      END DO
	    ELSE
	      IKT_ST=1
	      DO IKT=1,NKT
	        IF(XKT(IKT) .GE. THD(IT)%XTAB(1))THEN
	          IKT_ST=IKT
	          EXIT
	        END IF
	      END DO
	      IKT_END=NKT
	      DO IKT=IKT_ST,NKT
	        IF(XKT(IKT) .GE. THD(IT)%XTAB(THD(IT)%NTAB))THEN
	          IKT_END=IKT-1
	          EXIT
	        END IF
	      END DO
!
	      J=1
	      DO IKT=IKT_ST,IKT_END
	        DO WHILE(XKT(IKT) .GT. THD(IT)%XTAB(J+1))
	          J=J+1
	        END DO
	        T1=DLOG(XKT(IKT)/THD(IT)%XTAB(J))/DLOG(THD(IT)%XTAB(J+1)/THD(IT)%XTAB(J))
	        THD(IT)%CROSS_SEC(IKT)=EXP( T1*DLOG(THD(IT)%YTAB(J+1)) + (1.0D0-T1)*DLOG(THD(IT)%YTAB(J)) ) 
	      END DO
	      I=THD(IT)%NTAB
	      T1=DLOG(THD(IT)%YTAB(I)/THD(IT)%YTAB(I-1))/DLOG(THD(IT)%XTAB(I)/THD(IT)%XTAB(I-1))
	      DO IKT=IKT_END+1,NKT
	        THD(IT)%CROSS_SEC(IKT)=THD(IT)%YTAB(I)*(XKT(IKT)/THD(IT)%XTAB(I))**T1
	      END DO
	    END IF
	    THD(IT)%CROSS_SEC=THD(IT)%CROSS_SEC*ATM(ID)%ION_CROSEC_NTFAC
!
	    IF(NEGATIVE_CROSS_SEC)THEN
	      ID=THD(IT)%LNK_TO_ION
	      WRITE(6,*)'Negative cross section(s) for ',ION_ID(ID),' in ARNAUD_CROSS_V2!'
	      WRITE(6,*)'Forcing them to be zeros.'
	    END IF
	  END IF
	END DO
!	
	RETURN
	END
