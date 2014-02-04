C
C  Routine to read in a R,V,GRADV grid from file  
C    R(ND)= RP / 10**10cm and V in km/seg and SIGMA is dimensionless
C    sigma=R/V dv/dr -1.
C
	SUBROUTINE STARRAVE(R,V,SIGMA,ND,LU,RMAX,RP)
        IMPLICIT NONE
C
C Altered 26-May-1996 ERROR_LU installed.
C
	INTEGER ND,LU,I,NDOLD
	REAL*8 R(ND),V(ND),SIGMA(ND)
	REAL*8 RMAX,RP
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	LUER=ERROR_LU()                        
        OPEN(UNIT=LU,STATUS='OLD',FILE='RDINR',ACTION='READ')
        READ(LU,*) NDOLD
C Check relative values.
        IF(ND .NE. NDOLD)THEN
   	  WRITE(LUER,*)'Error-NDOLD and ND are not equal in RDINR'
	  WRITE(LUER,*)'NDOLD=',NDOLD,' ND=',ND
          STOP
        END IF
        DO 100 I=1,ND
          READ(LU,*,END=200)R(I),V(I),SIGMA(I)
 100    CONTINUE
 200    CONTINUE
        RP=R(ND)
        RMAX=R(1)
        CLOSE(UNIT=LU)
        RETURN
 
	END
