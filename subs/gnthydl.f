C
C Route to return the gaunt factor for individual L states of hydrogen.
C The gaunt factors (g) for all l sub-levels of level n are returned. 
C KSQ is the scaled energy of the ejected electron in Rydbergs (==K*K)
C (i.e E/z/z).
C
C Based on HYDPHOT_SUB which in turn was based on Brockelhurst (171,
C MNRAS, 153, 471-490).
C

	SUBROUTINE GNTHYDL(SIGMA,G,N,KSQ)
	IMPLICIT NONE
C
C Altered 24-May-1996 - RONE etc inserted.
C Created 29-Aug-1989 - Verified against tbales of Karzas and Latter.
C
	INTEGER N,I1,I2
	REAL*8 SIGMA(0:N-1),KSQ
	REAL*8 G(0:N-1,0:N+1)
C
	REAL*8 GAMRAT
C
C Local variables
C
	INTEGER L
	REAL*8 G0,NKSQ,MULT,K,ROOT3
C
	REAL*8, PARAMETER :: PI=3.1415926535897932384D0
	REAL*8, PARAMETER :: RHALF=0.5D0
	REAL*8, PARAMETER :: RONE=1.0D0
	REAL*8, PARAMETER :: RTWO=2.0D0
	REAL*8, PARAMETER :: RTHREE=3.0D0
	REAL*8, PARAMETER :: RFOUR=4.0D0
C
	ROOT3=SQRT(RTHREE)
C
C NB - arguments of GAMRAT are integers. G0 is identical to the
C      definition of Burgess apart from a factor of SQRT( N*SQRT(3)/16 ).
C      The first SQRT arrises since the G values are squared to obtain
C      the cross sections.
C
	I1=2*N-1
	I2=1
	G0=SQRT(RTWO*ROOT3*PI*N)*N/GAMRAT(I1,I2) *
	1       ( (RFOUR*N)**N )*EXP(-RTWO*N)
C
	K=SQRT(KSQ)
	NKSQ=RONE+N*N*KSQ
	IF(K .EQ. 0)THEN
	  G(N-1,N)=G0
	ELSE
C
C Have multiplied by NKSQ^2 (i.e multiplied cross section by NKSQ^4).
C
	  G(N-1,N)=G0 * EXP(  RTWO*N -RTWO/K*ATAN( DBLE(N*K) )  )
	1         /( NKSQ**N ) / SQRT(RONE-EXP(-RTWO*PI/K))
	END IF
	IF( N .EQ. 1 )GOTO 1000
C
	G(N-2,N-1)=(RTWO*N-RONE)*NKSQ*N*G(N-1,N)
	G(N-1,N-2)=RHALF*NKSQ*G(N-1,N)/N
	IF( N .EQ. 2 )GOTO 1000
C
	G(N-2,N-3)=(RTWO*N-RONE)*( RFOUR+ (N-RONE)*NKSQ )*G(N-1,N-2)
C
	DO L=N-1,2,-1
	  G(L-2,L-1)=( RFOUR*N*N - RFOUR*L*L +
	1              L*(RTWO*L-RONE)*NKSQ )*G(L-1,L) -
	1              RFOUR*N*N*(N+L)*(N-L)*( RONE+(L+1)*(L+1)*KSQ )
	1              *G(L,L+1)
	END DO
	IF( N .EQ. 3 )GOTO 1000
C
	DO L=N-2,2,-1
	  G(L-1,L-2)=( RFOUR*N*N - RFOUR*L*L +
	1              L*(RTWO*L+RONE)*NKSQ )*G(L,L-1) -
	1              RFOUR*N*N*( N*N - (L+1)*(L+1) ) *
	1              (RONE + L*L*KSQ )*G(L+1,L)
	END DO
C
1000	CONTINUE
C
C Convert from old G to little g.
C
C Case A - LPRIME = L + 1.
C
	MULT=RONE
	DO L=0,N-1
	  MULT=MULT*(RONE+(L+1)*(L+1)*KSQ)
	  I1=N+L
	  I2=N-L-1
	  G(L,L+1)=SQRT( GAMRAT(I1,I2)*MULT )* ( (RTWO*N)**(L-N) )
	1            *G(L,L+1)
	END DO
C
C Case B - LPRIME = L - 1.
C
	MULT=RONE
	DO L=1,N-1
	  MULT=MULT*(RONE+(L-1)*(L-1)*KSQ)
	  I1=N+L
	  I2=N-L-1
	  G(L,L-1)=SQRT( GAMRAT(I1,I2)*MULT )* ( (RTWO*N)**(L-N) )
	1            *G(L,L-1)
	END DO
C
C Can now comput photoionization cross-sections. Note that the constant
C is included in the definition of G0 and G(N-1,N).
C
	DO L=1,N-1
	  SIGMA(L)=( L*G(L,L-1)*G(L,L-1)
	1                  +(L+RONE)*G(L,L+1)*G(L,L+1) )/(RTWO*L+RONE)
	END DO
	SIGMA(0)=G(0,1)*G(0,1)
C
	RETURN
	END
