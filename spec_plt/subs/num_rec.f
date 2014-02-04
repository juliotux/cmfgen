
C---------------------------------------------------------------------------
C The following routines are taken from Numerical Recipes
C---------------------------------------------------------------------------

      SUBROUTINE CONVLV(DATA,N,RESPNS,M,ISIGN,ANS)
      
      implicit none
      include 'constants.inc'
      include 'parameters.inc'

      integer*4 n,m,isign,i,no2
C      PARAMETER(NMAX=131072)
      real*8 data(NMAX), respns(NMAX)
      complex*16 FFT(NMAX),ANS(NMAX)
      

      DO 11 I=1,(M-1)/2
        RESPNS(N+1-I)=RESPNS(M+1-I)
11    CONTINUE
      DO 12 I=(M+3)/2,N-(M-1)/2
        RESPNS(I)=0.0
12    CONTINUE

      CALL TWOFFT(DATA,RESPNS,FFT,ANS,N)
      NO2=N/2
      DO 13 I=1,NO2+1
        IF (ISIGN.EQ.1) THEN
          ANS(I)=FFT(I)*ANS(I)/dble(NO2)
        ELSE IF (ISIGN.EQ.-1) THEN
          IF (CDABS(ANS(I)) .EQ. zero) PAUSE 'DECONVOLVING AT A RESPONSE ZERO'
          ANS(I)=FFT(I)/ANS(I)/dble(NO2)
        ELSE
          WRITE(6,*)'NO MEANING FOR ISIGN in CNVLV'
	  STOP
        ENDIF
13    CONTINUE
      ANS(1)=CMPLX(REAL(ANS(1)),REAL(ANS(NO2+1)),kind(0d0))
      CALL REALFT(ANS,NO2,-1)
      RETURN
      END

      SUBROUTINE TWOFFT(DATA1,DATA2,FFT1,FFT2,N)

      implicit none
      include 'parameters.inc'
      include 'constants.inc'

      integer*4 n
      real*8 DATA1(NMAX),DATA2(NMAX)
      COMPLEX*16 FFT1(NMAX),FFT2(NMAX),H1,H2,C1,C2
      integer*4 j,n2,i

      write(*,*)'twofft'

      C1=CMPLX(0.5d0,0.0d0,kind(0d0))
      C2=CMPLX(0.0d0,-0.5d0,kind(0d0))
      DO 11 J=1,N
        FFT1(J)=CMPLX(DATA1(J),DATA2(J),kind(0d0))
11    CONTINUE

      CALL FOUR1(FFT1,N,1)
      FFT2(1)=CMPLX(AIMAG(FFT1(1)),zero,kind(0d0))
      FFT1(1)=CMPLX(REAL(FFT1(1)),zero,kind(0d0))
      N2=N+2
      DO 12 J=2,N/2+1
        H1=C1*(FFT1(J)+CONJG(FFT1(N2-J)))
        H2=C2*(FFT1(J)-CONJG(FFT1(N2-J)))
        FFT1(J)=H1
        FFT1(N2-J)=CONJG(H1)
        FFT2(J)=H2
        FFT2(N2-J)=CONJG(H2)
12    CONTINUE
      RETURN
      END

      SUBROUTINE REALFT(DATA,N,ISIGN)

      implicit none
      include 'constants.inc'
      include 'parameters.inc'

      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      real*8 DATA(NMAX)
      integer*4 n,isign,i,i1,i2,i3,i4,n2p3
      real*8 wis,wrs
      real*8 c1,c2,h1r,h1i,h2r,h2i

      write(*,*)'realft'

      THETA=6.28318530717959D0/2.0D0/DBLE(N)
      C1=0.5d0
      IF (ISIGN.EQ.1) THEN
        C2=-0.5d0
        CALL FOUR1(DATA,N,+1)
      ELSE
        C2=0.5d0
        THETA=-THETA
      ENDIF
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2d0
      WPI=DSIN(THETA)
      WR=1.0D0+WPR
      WI=WPI
      N2P3=2*N+3

      DO 11 I=2,N/2+1
        I1=2*I-1
        I2=I1+1
        I3=N2P3-I2
        I4=I3+1
        WRS=(WR)
C        WRS=SNGL(WR)
C        WIS=SNGL(WI)
        WIS=(WI)
        H1R=C1*(DATA(I1)+DATA(I3))
        H1I=C1*(DATA(I2)-DATA(I4))
        H2R=-C2*(DATA(I2)+DATA(I4))
        H2I=C2*(DATA(I1)-DATA(I3))
        DATA(I1)=H1R+WRS*H2R-WIS*H2I
        DATA(I2)=H1I+WRS*H2I+WIS*H2R
        DATA(I3)=H1R-WRS*H2R+WIS*H2I
        DATA(I4)=-H1I+WRS*H2I+WIS*H2R
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
11    CONTINUE
      IF (ISIGN.EQ.1) THEN
        H1R=DATA(1)
        DATA(1)=H1R+DATA(2)
        DATA(2)=H1R-DATA(2)
      ELSE
        H1R=DATA(1)
        DATA(1)=C1*(H1R+DATA(2))
        DATA(2)=C1*(H1R-DATA(2))
        CALL FOUR1(DATA,N,-1)
      ENDIF
      RETURN
      END

      SUBROUTINE FOUR1(DATA,NN,ISIGN)

      implicit none
      include 'constants.inc'
      include 'parameters.inc'

      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      integer*4 N,NN,I,J,M,mmax,istep,isign
      real*8 data(NMAX)
      real*8 tempi,tempr

      write(*,*)'four1'

      
      N=2*NN

      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/dble(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2d0
        WPI=DSIN(THETA)
        WR=one
        WI=zero
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=(WR)*DATA(J)-(WI)*DATA(J+1)
            TEMPI=(WR)*DATA(J+1)+(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END
