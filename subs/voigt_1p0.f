C
C***************************************************************
C
	FUNCTION VOIGT(A,V)
C
C Subroutine to evaluate a VOIGT profile, normalized so that the integral
C over V is unity.
C
C Routine is designed for speed.
C
C Adapted from a Midas which was based on a routined obtained from
C  S L Wright /UCL/ and a routine from Peter Hofflich which inturn
C was based on HUI,A.K.,ARMSTRONG,E.H.,WRAY,A.A. 1978,JQSRT 19,P.509
C
C Accuracy:
C         < 0.1% for 0.1 < a 
C         < 0.5% for a < 0.1
C         < 0.3% for a < 0.05
C
	IMPLICIT NONE
C                                                                   
	REAL*8  VOIGT,V,A
C
C Variables for interpolation (a < 0.1).
C
	INTEGER N,N1
	REAL*8  V0,V2
	REAL*8  H1(81)
C
C variables for functional form (a > 0.1)
C
	COMPLEX ZH,F
	REAL*8 C(7),D(7)
C
C Interpolation constants.
C
	DATA H1/
	1  -1.1283800D0,-1.1059600D0,-1.0404800D0,-0.9370300D0,-0.8034600D0,-0.6494500D0,
	1  -0.4855200D0,-0.3219200D0,-0.1677200D0,-0.0301200D0, 0.0859400D0, 0.1778900D0,
	1   0.2453700D0, 0.2898100D0, 0.3139400D0, 0.3213000D0, 0.3157300D0, 0.3009400D0,
	1   0.2802700D0, 0.2564800D0, 0.2317260D0, 0.2075280D0, 0.1848820D0, 0.1643410D0,
	1   0.1461280D0, 0.1302360D0, 0.1165150D0, 0.1047390D0, 0.0946530D0, 0.0860050D0,
	1   0.0785650D0, 0.0721290D0, 0.0665260D0, 0.0616150D0, 0.0572810D0, 0.0534300D0,
	1   0.0499880D0, 0.0468940D0, 0.0440980D0, 0.0415610D0, 0.0392500D0, 0.0351950D0,
	1   0.0317620D0, 0.0288240D0, 0.0262880D0, 0.0240810D0, 0.0221460D0, 0.0204410D0,
	1   0.0189290D0, 0.0175820D0, 0.0163750D0, 0.0152910D0, 0.0143120D0, 0.0134260D0,
	1   0.0126200D0, 0.0118860D0, 0.0112145D0, 0.0105990D0, 0.0100332D0, 0.0095119D0,
	1   0.0090306D0, 0.0085852D0, 0.0081722D0, 0.0077885D0, 0.0074314D0, 0.0070985D0,
	1   0.0067875D0, 0.0064967D0, 0.0062243D0, 0.0059688D0, 0.0057287D0, 0.0055030D0,
	1   0.0052903D0, 0.0050898D0, 0.0049006D0, 0.0047217D0, 0.0045526D0, 0.0043924D0,
	1   0.0042405D0, 0.0040964D0, 0.0039595D0/
C
C Functional constants.
C
	DATA C/122.6079D0,214.3823D0,181.9285D0,93.15558D0,30.18014D0,5.912626D0,5.641896D-01/
	DATA D/122.6079D0,352.7306D0,457.3344D0,348.7039D0,170.3540D0,53.99291D0,10.47986D0/
C
	IF(A .GT. 0.1D0)THEN
          ZH = CMPLX(A,(-ABS(V)))
          F = ((((((C(7)*ZH+C(6))*ZH+C(5))*ZH+C(4))*ZH+C(3))*ZH+C(2))
	1      *ZH+C(1))
	1      /(((((((ZH+D(7))*ZH+D(6))*ZH+D(5))*ZH+D(4))*ZH+D(3))*ZH+D(2)
	1      )*ZH+D(1))
	  VOIGT = REAL(F)
	ELSE IF(A .EQ. 0.0D0)THEN
	  VOIGT=EXP(-V*V)
	ELSE
C
C In original UCL routine, H0 and H2 were interpolated. Better to compute ---
C not much slower and improves accuracy from 3% to beter than 0.5% for
C a=0.001.
C
	  V0=ABS(V)*10.D0
	  N=V0
          IF(N.LT.40)THEN
	    V2=V0-N                     
	    N=N+1
	    N1=N+1
	    VOIGT=(1.0D0+A*A*(1.0D0-2*V*V))*EXP(-V*V) + 
	1            A*(H1(N)+V2*(H1(N1)-H1(N)) )
	  ELSE IF(N.LT.120)THEN
	    N=N/2+20
     	    V2=20.0D0+0.5D0*V0-N
	    N=N+1
	    N1=N+1
	    VOIGT=A*((H1(N1)-H1(N))*V2+H1(N))
	    IF(V0 .LT. 50.0D0)VOIGT=VOIGT+EXP(-V*V)
	  ELSE
            VOIGT=(0.56419D0+0.846D0/(V*V))/(V*V)*A
	  END IF
	END IF
	VOIGT=VOIGT/1.772453851D0	  	!SQRT(PI)
C
	RETURN
	END
