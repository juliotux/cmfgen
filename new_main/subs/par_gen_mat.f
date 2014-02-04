!
! Map the small BA rray onto the full BA array. The following works when
! the change in  /\z=1
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    DO J=1,SE(ID)%N_IV
	      JJ=SE(ID)%LNK_TO_F(J)
	      DO I=ATM(ID)%NXzV+1,SE(ID)%N_SE-1
	        C_ION(ID,JJ)=C_ION(ID,JJ)+SE(ID)%BA(I,J,BAND_INDX,K)
	      END DO
	    END DO
	  END DO
	END DO
!
! The following is when Xrays are included, and which change z by 2
!
	IF(XRAYS)THEN
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-2
              I=SE(ID)%XRAY_EQ
	      IF(I .NE. 0)THEN
	        DO J=1,SE(ID)%N_SE
	          JJ=SE(ID)%LNK_TO_F(J)
	          C_ION(ID+1,JJ)=C_ION(ID+1,JJ)+SE(ID)%BA(I,J,BAND_INDX,K)
	          C_ION(ID+2,JJ)=C_ION(ID+2,JJ)+SE(ID)%BA(I,J,BAND_INDX,K)
	        END DO
	      END IF
	    END DO
	  END DO
	END IF
!
! The following is completely general. No assumptions about how
! state changes ionization state.
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    DO 
	    DO J=1,SE(ID)%N_IV
	      JJ=SE(ID)%LNK_TO_F(J)
	      DO I=ATM(ID)%NXzV+1,SE(ID)%N_SE-1
	        C_ION(ID,JJ)=C_ION(ID,JJ)+SE(ID)%BA(I,J,BAND_INDX,K)
	        NID=SE(ID)%ION_EQ_ID(I)
	        C_ION(NID,JJ)=C_ION(NID,JJ)-SE(ID)%BA(I,J,BAND_INDX,K)
	      END DO
	    END DO
	  END DO
	END DO
!
