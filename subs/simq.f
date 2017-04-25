C        SUBROUTINE SIMQ                                               
C                                                                      
C        PURPOSE                                                       
C           OBTAIN SOLUTION OF A SET OF SIMULTANEOUS LINEAR EQUATIONS, 
C           AX=B                                                       
C                                                                      
C        USAGE                                                         
C           CALL SIMQ(A,B,N,KS)                                        
C                                                                      
C        DESCRIPTION OF PARAMETERS                                     
C           A - MATRIX OF COEFFICIENTS STORED COLUMNWISE.  THESE ARE   
C               DESTROYED IN THE COMPUTATION.  THE SIZE OF MATRIX A IS 
C               N BY N.                                                
C           B - VECTOR OF ORIGINAL CONSTANTS (LENGTH N). THESE ARE     
C               REPLACED BY FINAL SOLUTION VALUES, VECTOR X.           
C           N - NUMBER OF EQUATIONS AND VARIABLES. N MUST BE .GT. ONE. 
C           KS - OUTPUT DIGIT                                          
C                0 FOR A NORMAL SOLUTION                               
C                1 FOR A SINGULAR SET OF EQUATIONS                     
C                                                                      
C        REMARKS                                                       
C           MATRIX A MUST BE GENERAL.                                  
C           IF MATRIX IS SINGULAR , SOLUTION VALUES ARE MEANINGLESS.   
C           AN ALTERNATIVE SOLUTION MAY BE OBTAINED BY USING MATRIX    
C           INVERSION (MINV) AND MATRIX PRODUCT (GMPRD).               
C                                                                      
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                 
C           NONE                                                       
C                                                                      
C        METHOD                                                        
C           METHOD OF SOLUTION IS BY ELIMINATION USING LARGEST PIVOTAL 
C           DIVISOR. EACH STAGE OF ELIMINATION CONSISTS OF INTERCHANGIN
C           ROWS WHEN NECESSARY TO AVOID DIVISION BY ZERO OR SMALL     
C           ELEMENTS.                                                  
C           THE FORWARD SOLUTION TO OBTAIN VARIABLE N IS DONE IN       
C           N STAGES. THE BACK SOLUTION FOR THE OTHER VARIABLES IS     
C           CALCULATED BY SUCCESSIVE SUBSTITUTIONS. FINAL SOLUTION     
C           VALUES ARE DEVELOPED IN VECTOR B, WITH VARIABLE 1 IN B(1), 
C           VARIABLE 2 IN B(2),........, VARIABLE N IN B(N).           
C           IF NO PIVOT CAN BE FOUND EXCEEDING A TOLERANCE OF 0.0,     
C           THE MATRIX IS CONSIDERED SINGULAR AND KS IS SET TO 1. THIS 
C           TOLERANCE CAN BE MODIFIED BY REPLACING THE FIRST STATEMENT.
C                                                                      
C     .................................................................
C                                                                      
      SUBROUTINE SIMQ(A,B,N,KS)
      IMPLICIT NONE
C
C Altered 05-Dec-1996: Aritmetic if's  removed. END DO used to terminate all
C                         DO LOOPS. Changes verified.
C Altered 26-May-1996: Implicit none installed.
C
      INTEGER N,KS
      REAL*8 A(N),B(N)
C
C Local variables.
C
      REAL*8 TOL
      REAL*8 BIGA
      REAL*8 SAVE
C
      INTEGER I,J,K
      INTEGER I1,I2,IJ,IT,IA,IB,IC,IX
      INTEGER IQS,IMAX,IXJ,IXJX
      INTEGER JJ,JY,JX,JJX,NY
C
C Forward solution
C 
      TOL=0.0
      KS=0   
      JJ=-N  
      DO J=1,N
        JY=J+1
        JJ=JJ+N+1
        BIGA=0.0D0
        IT=JJ-J
        DO I=J,N
C
C Search for maximum coefficient in column
C                                                                      
          IJ=IT+I                                                          
          IF(ABS(BIGA) .LT. ABS(A(IJ)))THEN
            BIGA=A(IJ)
            IMAX=I
          END IF
        END DO
C                                      
C Test for pivot less than tolerance (singular matrix)          
C                                                                      
        IF(ABS(BIGA) .LE. TOL)THEN
          KS=1                                                             
          RETURN                                                           
        END IF
C
C Interchange rows if necessary.
C
        I1=J+N*(J-2)                                                     
        IT=IMAX-J                                                        
        DO K=J,N                                                      
          I1=I1+N                                                          
          I2=I1+IT                                                         
          SAVE=A(I1)                                                       
          A(I1)=A(I2)                                                      
          A(I2)=SAVE
C
C Divide equation by leading coefficient.
C 
          A(I1)=A(I1)/BIGA                                                 
        END DO
        SAVE=B(IMAX)                                                     
        B(IMAX)=B(J)                                                     
        B(J)=SAVE/BIGA                                                   
C                                                                      
C        ELIMINATE NEXT VARIABLE                                       
C                                                                      
        IF(J .NE. N)THEN
          IQS=N*(J-1)               
          DO IX=JY,N                   
            IXJ=IQS+IX                      
            IT=J-IX                         
            DO JX=JY,N                   
              IXJX=N*(JX-1)+IX                
              JJX=IXJX+IT                     
              A(IXJX)=A(IXJX)-(A(IXJ)*A(JJX)) 
	    END DO
            B(IX)=B(IX)-(B(J)*A(IXJ))       
	  END DO
        END IF
      END DO
C                                                             
C        BACK SOLUTION                                        
C                                                             
      NY=N-1        
      IT=N*N        
      DO J=1,NY  
        IA=IT-J       
        IB=N-J        
        IC=N          
        DO K=1,J   
          B(IB)=B(IB)-A(IA)*B(IC)  
          IA=IA-N                                                 
          IC=IC-1  
        END DO
      END DO
C
      RETURN   
      END      
