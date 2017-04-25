!===================================================================
!
      SUBROUTINE RECOM_OPAC_V2(SIGMA,NU,EDGE,FREQ_SCL_FAC,STAT_WT,GION,NPHOT,
     *        np_max,TOT_REC,TEMP)
!
!===================================================================
!
! Routine to compute the recombination coefficients using
! photionization cross-sections from the opacity project. Stimulated
! recombination is NOT taken into account. The fequency should
! be in units of the the EDGE frequency. EDGE frequency should be in
! units of 10^15 Hz. A combination of Simpons rule, and the trapazoidal
! ruke is used.
!
! Created as separate subroutine 10-Oct-1990 : GION installed
!
! Altered 10-Oct-1990 : Intgeration section changed. Now checks whether
!                       sufficient points are available to accurately
!                       integrate across exponential. If not, additional
!                       points are inserted. Old routine was giving
!                       inaccurate answers for opacity grid, especially
!                       above 2p threshold (eg 2p3d3Fo state).
!
! Altered   22/6/96 DLM Added parameter scratch so wouldn't have to pass
!                       variables cross and rec.
!                       Also installed check so integration stops if
!                       the exponent of exp(-h*nu/k/t) is greater than
!                       -100.  This solves the problem of getting a floating
!                       point operation error when this exponentation becomes
!                       very small.  The outputted recombination coefficients
!                       are not used anyway.
!
! Altered  16-Jun-1999 DJH: Routine was giving in accurate answers when
!                        NU(IST) was note exactly 1. Needed to allow for
!                        the insertion of extra points. Introduced variables
!                        LST_CROSS (replaces CROSS(JEND) and LST_NU so that
!                        the same ste of control statements could be used for
!                        all frequency values.
!
! Altered   2-Nov-2009 DJH: Call changed. 
!                           Adjusted to have alternative frequenct scaling.
!                           Done because of photoiozation of states above ionization limit.
!
	IMPLICIT NONE
!
        integer scratch
        parameter(scratch=100000)
!
	INTEGER NPHOT,np_max
	REAL*8 SIGMA(np_max),NU(np_max),STAT_WT,GION,TEMP
!
	REAL*8 REC(scratch),CROSS(scratch),TOT_REC
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	INTEGER I,J,NJ,IST
	REAL*8 T1,TP1,CONST,T2
	REAL*8 NU_STEP,DEL_NU,DEL_NU_2,DEL_REC
	REAL*8 FRAC_DIFF,MIN_SPACING,X,SIG
	REAL*8 EDGE,SIG_THRESH
	REAL*8 FREQ_SCL_FAC
	REAL*8 LST_NU
	REAL*8 LST_CROSS
	LOGICAL EQUAL
	EXTERNAL EQUAL
!
	CHIBF=2.815E-06
	CHIFF=3.69E-29
	HDKT=4.7994145
	TWOHCSQ=0.0147452575
	OPLIN=2.6540081E+08
	EMLIN=5.27296E-03
!
        if(scratch.lt.np_max)then
          print*,' parameter scratch .lt. np_max'
          print*,' scratch =',scratch
          print*,'  np_max =',np_max
          stop ' stopped in `recom_opac`'
        endif
!
! FRAC_DIFF is the accuracy to which adjacent step sizes must be equal
! is SImpson rule is to be used.
!
	FRAC_DIFF=1.0D-06
!
! MIN_SPACING determines the maximum spacing, before extra points are
! inserted to resolve the exponential function. Defined so that 0.2
! is the maximum step in the exponent of the exponential between adjacent
! grid points. (i.e. exp(0), exp(-0.2), exp(-0.4)  etc.).
!
	MIN_SPACING=0.2*TEMP/HDKT/FREQ_SCL_FAC		!EDGE
!
        if(nphot.gt.np_max)then
          print*,' increase np_max to .gt. nphot,',nphot
          stop ' stopped in recom_opacity'
        endif
!
	TOT_REC=0.0D0
	DO I=1,NPHOT
	  REC(I)=0.0D0
	END DO
!
! Determine index of first frequency above (or equal to) threshold.
! Determine cross-section at threshold.
!
	IST=1
	DO I=1,NPHOT
	  IF(NU(I) .LT. 1.0D0)IST=IST+1
	END DO
	IF(IST .NE. 1)THEN
	  T1=(NU(IST)-1.0D0)/(NU(IST)-NU(IST-1))
	  SIG_THRESH=(1.0D0-T1)*SIGMA(IST)+T1*SIGMA(IST-1)
	ELSE
	  SIG_THRESH=SIGMA(1)
	END IF
!
	TP1=1.8965D+17		!4*PI/H*DEX(-10)
	TP1=TP1*1.0D-08         !As sigma in megabarns.
!
! The factor EDGE**3 appears in the CONST definition as NU is in units
! of EDGE (**2 from nu^2 and additional factor from dnu).
!
	T2=HDKT*FREQ_SCL_FAC/TEMP
	CONST=TP1*TWOHCSQ*(FREQ_SCL_FAC**3)*2.07D-22*STAT_WT/GION
	CONST=CONST*EXP(HDKT*(EDGE-FREQ_SCL_FAC)/TEMP)
!
! We can now compute recombination rate. We first allow for the
! recombination near threshold.
!
! We use lST_NU and LST_CROSS to access the integrand in the CROSS array
! at the end point of the previous integrations. Saves some computational
! effort.
!
! NB: NU=1 at threshold, and thus [NU^2 Exponential term] is unity.
!
	LST_CROSS=SIG_THRESH*CONST
	LST_NU=1.0D0
	I=IST
	IF(IST .GT. 1)I=I-1		!To get bit between 1 and NU(IST)
	DO WHILE (I .LT. NPHOT)
!
! If exponent is -500 then can quit calculating recombination coefficient.
! This was installed to fix "floating point operation error" caused by taking the
! exponatial of a number of order -500.
!
          if(-t2*(nu(i+1)-1.0D0) .lt.-500.)then
            rec(nphot)=rec(i)
            goto 100
          endif
!
	  DEL_NU=NU(I+1)-LST_NU
	  IF( DEL_NU .LT. MIN_SPACING)THEN
	    IF(I .LT. NPHOT-1)THEN
	      DEL_NU_2=NU(I+2)-NU(I+1)
	    ELSE
	      DEL_NU_2=0.0D0
	    END IF
	    IF( EQUAL(DEL_NU,DEL_NU_2,FRAC_DIFF) )THEN
!
! Use Simpsons rule.
!
	      CROSS(1)=LST_CROSS
	      CROSS(2)=CONST*( NU(I+1)**2 )*
     *             EXP( -T2*(NU(I+1)-1.0D0) )*SIGMA(I+1)
	      CROSS(3)=CONST*( NU(I+2)**2 )*
     *             EXP( -T2*(NU(I+2)-1.0D0) )*SIGMA(I+2)
	      DEL_REC=DEL_NU*(CROSS(1)+4.0D0*CROSS(2)+CROSS(3))/3.0D0
	      REC(I+2)=REC(I)+DEL_REC
	      REC(I+1)=REC(I)+0.5D0*DEL_REC		!Estimate only.
	      LST_CROSS=CROSS(3)
	      LST_NU=NU(I+2)
	      I=I+2
	    ELSE
!
! Use Trapazoidal rule.
!
	      CROSS(1)=LST_CROSS
	      CROSS(2)=CONST*( NU(I+1)**2 )*
     *             EXP( -T2*(NU(I+1)-1.0D0) )*SIGMA(I+1)
	      DEL_REC=DEL_NU*(CROSS(1) + CROSS(2))/2.0D0
	      REC(I+1)=REC(I)+DEL_REC
	      LST_CROSS=CROSS(2)
	      LST_NU=NU(I+1)
	      I=I+1
	    END IF
	  ELSE
!
! We need to insert extra points so that we can accurately integrate
! exponential weighting factor. We insert an odd number of points so that
! can use Simpsons rule.
!
	    NJ=INT(DEL_NU/MIN_SPACING)
	    IF( MOD(NJ,2) .EQ. 0)NJ=NJ+1
!
            if(nj.ge.np_max-1)then
              nj=np_max-2
              IF( MOD(NJ,2) .EQ. 0)NJ=NJ-1
            endif
!
	    NU_STEP=DEL_NU/(NJ+1)
	    CROSS(1)=LST_CROSS
!
! Point NJ+1 corresponds to NU(I+1).
!
	    DO J=1,NJ+1
	      X=LST_NU+NU_STEP*J
	      T1=J*NU_STEP/DEL_NU
	      SIG=T1*SIGMA(I+1)+(1.0D0-T1)*SIGMA(I)
	      CROSS(J+1)=CONST*( X**2 )*EXP( -T2*(X-1.0D0) )*SIG
	    END DO
!
! Can now perform integration.
!
	    DEL_REC=0
	    DO J=1,NJ,2
	      DEL_REC=DEL_REC +
     *            NU_STEP*(CROSS(J)+4.0D0*CROSS(J+1)+CROSS(J+2))/3.0D0
	    END DO
	    REC(I+1)=REC(I)+DEL_REC
	    LST_CROSS=CROSS(NJ+2)
	    LST_NU=NU(I+1)
	    I=I+1
	  END IF
	END DO
!
 100    TOT_REC=REC(NPHOT)
	IF(TOT_REC .NE. 0.0D0)THEN
	  DO I=IST,NPHOT
	    REC(I)=REC(I)/TOT_REC
	    WRITE(100,*)I,NU(I),REC(I)
	  END DO
	END IF
!
! Correct TOT_REC for frequency independent temperature dependance.
!
	TOT_REC=TOT_REC/(TEMP**1.5)
!
	RETURN
	END
