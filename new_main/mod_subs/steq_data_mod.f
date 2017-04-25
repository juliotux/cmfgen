!
! Module containing basic data for:
!                   (i) Each species
!          	   (ii) Model atom data (populations etc)
!		  (iii) Basic atmospheric structure.
!
      MODULE STEQ_DATA_MOD
!
!Altered 19-Aug-2015: BA_MAX_IONS_PER_SPECIES increaed to 21 (cur_hmi, 12-Jun-2015)
!                            
! Number of atomic species (e.g. H, C, N is 3 species).
!
	INTEGER, PARAMETER :: BA_NUM_SPECIES=26
!
! Maximum number of ionization stages per species. For H, need at this number
! has to be 2 or higher (as I and II). A setting of 10 implies that we can treat
! full atoms for ION_IX.
!
	INTEGER, PARAMETER :: BA_MAX_IONS_PER_SPECIES=21
	INTEGER, PARAMETER :: BA_MAX_NUM_IONS=BA_NUM_SPECIES*BA_MAX_IONS_PER_SPECIES
!
! Maximum number of photoionization routes for each species.
!
	INTEGER, PARAMETER :: BA_NPHOT_MAX=4
!
! Storage locations for the Statistical euilibrium equations, and its
! variation.
!
	TYPE STAT_EQ_DATA 
          REAL*8, POINTER :: STEQ(:,:)  		!Statistical equilibrium Eqns. for XzV
          REAL*8, POINTER :: STEQ_ADV(:)  		!Statistical equilibrium ion eqns. for advection terms
          REAL*8, POINTER :: QFV_P(:,:)  		!
          REAL*8, POINTER :: QFV_R(:,:)  		!
          REAL*8, POINTER :: BA(:,:,:,:)		!BA matrix for XzV
          REAL*8, POINTER :: BA_PAR(:,:,:)		!BA matrix for XzV (diagonal terms)
	  INTEGER, POINTER :: LNK_TO_IV(:)    	!
	  INTEGER, POINTER :: LNK_TO_F(:)     	!
	  INTEGER, POINTER :: EQ_IN_BA(:)     	!
	  INTEGER, POINTER :: EQ_TO_ION_LEV_PNT(:)	!
	  INTEGER, POINTER :: ION_LEV_TO_EQ_PNT(:)	!
	  INTEGER, POINTER :: STRT_ADV_ID(:)
	  INTEGER, POINTER :: END_ADV_ID(:)
	  INTEGER N_SE                        	!Number of S.E. Eqns. for XzV
	  INTEGER N_IV                        	!Number of important variables for XzV
	  INTEGER NUMBER_BAL_EQ			!Conservation equation
	  INTEGER XRAY_EQ				!Equation for Auger ionizations
	  LOGICAL Xzv_PRES                      	!Indicates whether ion is present.
          LOGICAL IMPURITY_SPECIES
	END TYPE STAT_EQ_DATA
!
        REAL*8, ALLOCATABLE :: STEQ_T(:)
        REAL*8, ALLOCATABLE :: STEQ_ED(:)
        REAL*8, ALLOCATABLE :: BA_T(:,:,:)
        REAL*8, ALLOCATABLE :: BA_ED(:,:,:)
        REAL*8, ALLOCATABLE :: BA_ADV_TERM(:,:)
!
        REAL*8, ALLOCATABLE :: BA_T_PAR(:,:)
!
        TYPE (STAT_EQ_DATA)  SE(BA_NUM_SPECIES*BA_MAX_IONS_PER_SPECIES)
!
      END MODULE STEQ_DATA_MOD
