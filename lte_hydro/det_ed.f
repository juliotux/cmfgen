	SUBROUTINE DET_ED(ND,ABUND_SUM,DO_LEV_DIS)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
	INTEGER ND
	REAL*8 ABUND_SUM
	LOGICAL DO_LEV_DIS
!
	REAL*8 ION_POPS(ND,NUM_IONS)
	REAL*8 FSAHA(NUM_IONS)
	REAL*8 XZ(NUM_IONS)
	REAL*8 XZW(NUM_IONS)
	REAL*8 ED_OLD(ND)
	REAL*8 POPION_OLD(ND)
!
	REAL*8 XEDW
	REAL*8 XED_OLD
	REAL*8 XERR
	REAL*8 XKBT
	REAL*8 XG0
	REAL*8 XG1
	REAL*8 XGE
	REAL*8 T1,T2,T3
	REAL*8 PI
	REAL*8 HSQR
	REAL*8, PARAMETER :: HDKT=4.7994145D0
	REAL*8, PARAMETER :: EPS=1.0D-05
!
	REAL*8 PLANCKS_CONSTANT
	REAL*8 BOLTZMANN_CONSTANT
	REAL*8 ELECTRON_MASS
	EXTERNAL PLANCKS_CONSTANT 
	EXTERNAL BOLTZMANN_CONSTANT
	EXTERNAL ELECTRON_MASS
!
	INTEGER ISTART,IEND
	INTEGER L
	INTEGER K
	INTEGER I
	INTEGER IW
	INTEGER ID
	INTEGER ISPEC
	LOGICAL CONVERGED
!
	PI=ACOS(-1.0D0)
        HSQR = PLANCKS_CONSTANT()*PLANCKS_CONSTANT()
	DO ISPEC=1,NUM_SPECIES
	  DO L=1,ND
            POP_SPECIES(L,ISPEC)=AT_ABUND(ISPEC)*POP_ATOM(L)/ABUND_SUM ! 1/cm^3
          END DO
	END DO
!
! Perform initializations.
!
	ED_OLD=0.0D0
	ED=POP_ATOM
	POPION=POP_ATOM
	POPION_OLD=0.0D0
	FSAHA=0.0D0
	CONVERGED=.FALSE.
!
! Compute the effective statistical weight for all levels in all ions.
! The effective statistical weight is the product of the actual statistical
! weight, and the Boltzmann excitation factor. We store it in the FULL
! level population vector.
!
	DO L=1,ND
          XKBT = BOLTZMANN_CONSTANT() * T(L) * 1.0D4 
          T1 = (2.0D0*PI*ELECTRON_MASS()*XKBT/HSQR)**1.5
          T2 = HDKT/T(L)
          DO ID=1,NUM_IONS
            IF (ATM(ID)%XzV_PRES) THEN
              ATM(ID)%XzV_F(1,L)= ATM(ID)%GXzV_F(1)
	      DO I=2,ATM(ID)%NXzV_F
	        ATM(ID)%XzV_F(I,L)=ATM(ID)%GXzV_F(I)*EXP(T2*(ATM(ID)%EDGEXzV_F(I)-ATM(ID)%EDGEXzV_F(1)))
	      END DO
            END IF
          END DO
	END DO
!
	DO WHILE(.NOT. CONVERGED)
	   CONVERGED=.TRUE.
!
! Compute the occupation probabilities, first updating the dissolution constants.
!
	  CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DIS,ND)
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      CALL OCCUPATION_PROB(ATM(ID)%W_XZV_F,ATM(ID)%EDGEXzV_F,
	1       ATM(ID)%ZXzV,ATM(ID)%NXzV_F,ND)
	    END IF
	  END DO
!
! Solve the Saha equation for the mixture considered: we obtain the electron density
! and the relative population of all ions.
!
	  DO L=1,ND
	    DO WHILE( ABS(1.0D0-ED_OLD(L)/ED(L)) .GT. EPS .OR.  
	1             ABS(1.0D0-POPION_OLD(L)/POPION(L)) .GT. EPS )
	      ED_OLD(L)=ED(L)
	      POPION_OLD(L)=POPION(L)
	      CONVERGED=.FALSE.
!
              XKBT = BOLTZMANN_CONSTANT() * T(L) * 1.0D4 
              T1 = (2.0D0*PI*ELECTRON_MASS()*XKBT/HSQR)**1.5
              T2 = HDKT/T(L)
              DO ID=1,NUM_IONS
                IF (ATM(ID)%XzV_PRES) THEN
                  XGE = 2.0D0
                  XG0 = ATM(ID)%GXzV_F(1)
	          DO I=2,ATM(ID)%NXzV_F
	            XG0=XG0+ATM(ID)%W_XzV_F(I,L)*ATM(ID)%XzV_F(I,L)
	          END DO
                  XG1 = ATM(ID)%GIONXzV_F
                  IF(ATM(ID+1)%XzV_PRES) THEN
	            DO I=2,ATM(ID+1)%NXzV_F
	              XG1=XG1+ATM(ID+1)%W_XzV_F(I,L)*ATM(ID+1)%XzV_F(I,L)
                    END DO
	          END IF
	          FSAHA(ID) = (XG1*XGE/XG0) * T1 * DEXP(-ATM(ID)%EDGEXzV_F(1)*T2)  
                END IF
              END DO	!ion loop.
!
	      XEDW = ED(L)
              XED_OLD = XEDW
              XERR = 2.0D0 * EPS
	      DO WHILE (XERR .GT. EPS)
                T3 = 0.0D0                
                DO ISPEC=1,NUM_SPECIES
                  ISTART = SPECIES_BEG_ID(ISPEC)
                  IF (ISTART.NE.0) THEN
                    IEND = SPECIES_END_ID(ISPEC)          
                    T1 = 1.0D0
                    T2 = 1.0D0 
                    DO IW=ISTART+1,IEND
                      T1 = (FSAHA(IW-1)/XEDW) * T1
                      T2 = T2 + T1
                      XZW(IW) = T1
                    END DO
                    XZ(ISTART) = 1.0D0 / T2
                    DO IW=ISTART+1,IEND
                      XZ(IW) = XZW(IW) * XZ(ISTART)
                      T3 = T3 + XZ(IW) * POP_SPECIES(L,ISPEC) * ATM(IW-1)%ZXzV
                    END DO
                    XZ(ISTART:IEND) = XZ(ISTART:IEND) * POP_SPECIES(L,ISPEC)
                  END IF
                END DO
                XED_OLD = XEDW
                XEDW = T3 
                XERR = DABS(1.0D0-XED_OLD/XEDW)
              END DO
              ED(L)=XEDW 
              ION_POPS(L,:) = XZ(:)  ! each entry corresponds to ions ID=1,NUM_IONS
!
	    END DO	!Only do depths that are not finished
	  END DO	!Loop over L index
!
! Get total ion population for level dissolution. Note that ATM(ID)$ZxzV
! is the charge on the ion AFTER the valence electron is removed
! (i.e., it =1 for HI).
!
	  POPION=0.0D0
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      IF(ATM(ID)%ZXzV .GT. 1.01D0)THEN
	        DO L=1,ND
	          POPION(L)=POPION(L)+ION_POPS(L,ID)
	        END DO
	      END IF
	    END DO
	    IF(SPECIES_BEG_ID(ISPEC) .NE. 0)THEN
	      DO L=1,ND
	        POPION(L)=POPION(L)+ION_POPS(L,SPECIES_END_ID(ISPEC))
	      END DO
	    END IF
	  END DO
!
	END DO		!All Ne are not accurate
!
! Compute the occupation probabilities, first updating the dissolution constants
!
	CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DIS,ND)
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    CALL OCCUPATION_PROB(ATM(ID)%W_XZV_F,ATM(ID)%EDGEXzV_F,
	1       ATM(ID)%ZXzV,ATM(ID)%NXzV_F,ND)
	  END IF
	END DO
!
! Compute the LTE populations of the levels in the full atom, taking level
! dissolution into account. We loop backwards over ID in order to have the
! ion population available.
!
	DO ID=NUM_IONS,1,-1
	  IF(ATM(ID)%XzV_PRES)THEN
!
	    DO L=1,ND
	      T1 = 0.0D0    
	      DO I=1,ATM(ID)%NXzV_F
	        T3 = HDKT*(ATM(ID)%EDGEXzV_F(1)-ATM(ID)%EDGEXzV_F(I))/T(L)
                T2 = ATM(ID)%W_XzV_F(I,L)*ATM(ID)%GXzV_F(I)*EXP(-T3)
                ATM(ID)%XzVLTE_F(I,L) = T2
                T1 = T1 + T2
              END DO
	      DO I=1,ATM(ID)%NXzV_F
	        ATM(ID)%XzVLTE_F(I,L) = ATM(ID)%XzVLTE_F(I,L) * ION_POPS(L,ID) / T1
	      END DO
	      IF(ATM(ID+1)%XzV_PRES)THEN
	         ATM(ID)%DXzV_F(L)=ATM(ID+1)%XzVLTE_F(1,L)
	      ELSE
	         ATM(ID)%DXzV_F(L)=ION_POPS(L,ID+1)
	      END IF
!
	      ATM(ID)%XzVLTE(:,L) = 0.0D0
	      DO I=1,ATM(ID)%NXzV_F
	        K=ATM(ID)%F_TO_S_XzV(I)
	        ATM(ID)%XzVLTE(K,L)=ATM(ID)%XzVLTE(K,L)+ATM(ID)%XzVLTE_F(I,L)
	      END DO
	      IF(ATM(ID+1)%XzV_PRES)THEN
	        ATM(ID)%DXzV(L)=ATM(ID+1)%XzVLTE(1,L)
	      ELSE
	        ATM(ID)%DXzV(L)=ION_POPS(L,ID+1)
	      END IF
	    END DO 		!Over depth
	  END IF   		!Ion present
	END DO			!Over ion
!
	RETURN
	END
