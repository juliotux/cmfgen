c===================================================================
c
      SUBROUTINE RECOM_OPAC(SIGMA,NU,EDGE,STAT_WT,GION,NPHOT,
     *        np_max,TOT_REC,TEMP)
c
c===================================================================
C
C Routine to compute the recombination coefficients using
C photionization cross-sections from the opacity project. Stimulated
C recombination is NOT taken into account. The fequency should
C be in units of the the EDGE frequency. EDGE frequency should be in
C units of 10^15 Hz. A combination of Simpons rule, and the trapazoidal
C ruke is used.
C
C Created as separate subroutine 10-Oct-1990 : GION installed
C
C Altered 10-Oct-1990 : Intgeration section changed. Now checks whether
C                       sufficient points are available to accurately
C                       integrate across exponential. If not, additional
C                       points are inserted. Old routine was giving
C                       inaccurate answers for opacity grid, especially
C                       above 2p threshold (eg 2p3d3Fo state).
C
c Altered   22/6/96 DLM Added parameter scratch so wouldn't have to pass
c                       variables cross and rec.
c                       Also installed check so integration stops if
c                       the exponent of exp(-h*nu/k/t) is greater than
c                       -100.  This solves the problem of getting a floating
c                       point operation error when this exponentation becomes
c                       very small.  The outputted recombination coefficients
c                       are not used anyway.
c
C Altered  16-Jun-1999 DJH: Routine was giving in accurate answers when
C                        NU(IST) was note exactly 1. Needed to allow for
C                        the insertion of extra points. Introduced variables
C                        LST_CROSS (replaces CROSS(JEND) and LST_NU so that
C                        the same ste of control statements could be used for
C                        all frequency values.
c
	IMPLICIT NONE
C
        integer scratch
        parameter(scratch=100000)
c
	INTEGER NPHOT,np_max
	REAL*8 SIGMA(np_max),NU(np_max),STAT_WT,GION,TEMP
c
	REAL*8 REC(scratch),CROSS(scratch),TOT_REC
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
C
	INTEGER I,J,NJ,IST
	REAL*8 T1,TP1,CONST,T2
	REAL*8 NU_STEP,DEL_NU,DEL_NU_2,DEL_REC
	REAL*8 FRAC_DIFF,MIN_SPACING,X,SIG
	REAL*8 EDGE,SIG_THRESH
	REAL*8 LST_NU
	REAL*8 LST_CROSS
	LOGICAL EQUAL
	EXTERNAL EQUAL
c
	CHIBF=2.815E-06
	CHIFF=3.69E-29
	HDKT=4.7994145
	TWOHCSQ=0.0147452575
	OPLIN=2.6540081E+08
	EMLIN=5.27296E-03
c
        if(scratch.lt.np_max)then
          print*,' parameter scratch .lt. np_max'
          print*,' scratch =',scratch
          print*,'  np_max =',np_max
          stop ' stopped in `recom_opac`'
        endif
C
C FRAC_DIFF is the accuracy to which adjacent step sizes must be equal
C is SImpson rule is to be used.
C
	FRAC_DIFF=1.0D-06
C
C MIN_SPACING determines the maximum spacing, before extra points are
C inserted to resolve the exponential function. Defined so that 0.2
C is the maximum step in the exponent of the exponential between adjacent
C grid points. (i.e. exp(0), exp(-0.2), exp(-0.4)  etc.).
C
	MIN_SPACING=0.2*TEMP/HDKT/EDGE
C
        if(nphot.gt.np_max)then
          print*,' increase np_max to .gt. nphot,',nphot
          stop ' stopped in recom_opacity'
        endif
c
	DO I=1,NPHOT
	  REC(I)=0.0D0
	END DO
C
C Determine index of first frequency above (or equal to) threshold.
C Determine cross-section at threshold.
C
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
C
	TP1=1.8965D+17		!4*PI/H*DEX(-10)
	TP1=TP1*1.0D-08         !As sigma in megabarns.
C
C The factor EDGE**3 appears in the CONST definition as NU is in units
C of EDGE (**2 from nu^2 and additional factor from dnu).
C
	T2=HDKT*EDGE/TEMP
	CONST=TP1*TWOHCSQ*(EDGE**3)*2.07D-22*STAT_WT/GION
C
C We can now compute recombination rate. We first allow for the
C recombination near threshold.
C
C We use lST_NU and LST_CROSS to access the integrand in the CROSS array
C at the end point of the previous integrations. Saves some computational
C effort.
C
C NB: NU=1 at threshold, and thus [NU^2 Exponential term] is unity.
C
	LST_CROSS=SIG_THRESH*CONST
	LST_NU=1.0D0
	I=IST
	IF(IST .GT. 1)I=I-1		!To get bit between 1 and NU(IST)
	DO WHILE (I .LT. NPHOT)
c
c If exponent is -500 then can quit calculating recombination coefficient.
c This was installed to fix "floating point operation error" caused by taking the
c exponatial of a number of order -500.
c
          if(-t2*(nu(i+1)-1.0D0) .lt.-500.)then
            rec(nphot)=rec(i)
            goto 100
          endif
c
	  DEL_NU=NU(I+1)-LST_NU
	  IF( DEL_NU .LT. MIN_SPACING)THEN
	    IF(I .LT. NPHOT-1)THEN
	      DEL_NU_2=NU(I+2)-NU(I+1)
	    ELSE
	      DEL_NU_2=0.0D0
	    END IF
	    IF( EQUAL(DEL_NU,DEL_NU_2,FRAC_DIFF) )THEN
C
C Use Simpsons rule.
C
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
C
C Use Trapazoidal rule.
C
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
C
C We need to insert extra points so that we can accurately integrate
C exponential weighting factor. We insert an odd number of points so that
C can use Simpsons rule.
C
	    NJ=INT(DEL_NU/MIN_SPACING)
	    IF( MOD(NJ,2) .EQ. 0)NJ=NJ+1
c
            if(nj.ge.np_max-1)then
              nj=np_max-2
              IF( MOD(NJ,2) .EQ. 0)NJ=NJ-1
            endif
c
	    NU_STEP=DEL_NU/(NJ+1)
	    CROSS(1)=LST_CROSS
C
C Point NJ+1 corresponds to NU(I+1).
C
	    DO J=1,NJ+1
	      X=LST_NU+NU_STEP*J
	      T1=J*NU_STEP/DEL_NU
	      SIG=T1*SIGMA(I+1)+(1.0D0-T1)*SIGMA(I)
	      CROSS(J+1)=CONST*( X**2 )*EXP( -T2*(X-1.0D0) )*SIG
	    END DO
C
C Can now perform integration.
C
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
C
 100    TOT_REC=REC(NPHOT)
	IF(TOT_REC .NE. 0.0D0)THEN
	  DO I=IST,NPHOT
	    REC(I)=REC(I)/TOT_REC
	  END DO
	END IF
C
C Correct TOT_REC for frequency independent temperature dependance.
C
	TOT_REC=TOT_REC/(TEMP**1.5)
C
	RETURN
	END
