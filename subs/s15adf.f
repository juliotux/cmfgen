C-----------------------------------------------------------------------
C   Purpose             - EVALUATE THE COMPLEMENTED ERROR FUNCTION OF   
C                           A DOUBLE PRECISION ARGUMENT                 
C                                                                       
C   Usage               - RESULT = S15ADF(Y)                             
C                                                                       
C   Arguments    Y      - Input double precision argument of the        
C                           complemented error function.                
C                ERFC  - Output double precision value of the          
C                           complemented error function. erfc must be  
C                           typed double precision in the calling       
C                           program.                                    
C
C-----------------------------------------------------------------------
C                                                                       
      FUNCTION S15ADF(Y,IFAIL)
C
C Specifications for arguments
C
      REAL*8 S15ADF,Y
C
C Local variables
C
      REAL*8 P(5),Q(4),P1(9),Q1(8),P2(6),Q2(5)              
      REAL*8 XMIN,XLARGE,SQRPI,X,           
     *                   RES,XSQ,XNUM,XDEN,XI,XBIG                      
      INTEGER            ISW,I,IFAIL                                          
C
C Coefficients for 0.0 .LE. Y .LT.     
C                                  .477                                 
      DATA               P(1)/113.8641541510502D0/,                     
     *                   P(2)/377.4852376853020D0/,                     
     *                   P(3)/3209.377589138469D0/,                     
     *                   P(4)/.1857777061846032D0/,                     
     *                   P(5)/3.161123743870566D0/                      
      DATA               Q(1)/244.0246379344442D0/,                     
     *                   Q(2)/1282.616526077372D0/,                     
     *                   Q(3)/2844.236833439171D0/,                     
     *                   Q(4)/23.60129095234412D0/                      
C                                  COEFFICIENTS FOR .477 .LE. Y         
C                                  .LE. 4.0                             
      DATA               P1(1)/8.883149794388376D0/,                    
     *                   P1(2)/66.11919063714163D0/,                    
     *                   P1(3)/298.6351381974001D0/,                    
     *                   P1(4)/881.9522212417691D0/,                    
     *                   P1(5)/1712.047612634071D0/,                    
     *                   P1(6)/2051.078377826071D0/,                    
     *                   P1(7)/1230.339354797997D0/,                    
     *                   P1(8)/2.153115354744038D-8/,                   
     *                   P1(9)/.5641884969886701D0/                     
      DATA               Q1(1)/117.6939508913125D0/,                    
     *                   Q1(2)/537.1811018620099D0/,                    
     *                   Q1(3)/1621.389574566690D0/,                    
     *                   Q1(4)/3290.799235733460D0/,                    
     *                   Q1(5)/4362.619090143247D0/,                    
     *                   Q1(6)/3439.367674143722D0/,                    
     *                   Q1(7)/1230.339354803749D0/,                    
     *                   Q1(8)/15.74492611070983D0/                     
C                                  COEFFICIENTS FOR 4.0 .LT. Y          
      DATA               P2(1)/-3.603448999498044D-01/,                 
     *                   P2(2)/-1.257817261112292D-01/,                 
     *                   P2(3)/-1.608378514874228D-02/,                 
     *                   P2(4)/-6.587491615298378D-04/,                 
     *                   P2(5)/-1.631538713730210D-02/,                 
     *                   P2(6)/-3.053266349612323D-01/                  
      DATA               Q2(1)/1.872952849923460D0/,                    
     *                   Q2(2)/5.279051029514284D-01/,                  
     *                   Q2(3)/6.051834131244132D-02/,                  
     *                   Q2(4)/2.335204976268692D-03/,                  
     *                   Q2(5)/2.568520192289822D0/                     
C                                  CONSTANTS                            
      DATA               XMIN/1.0D-10/,XLARGE/6.375D0/                  
C                                  ERFC(XBIG) .APPROX. DETAP           
      DATA               XBIG/13.3D0/                                   
      DATA               SQRPI/.5641895835477563D0/                     
C                                  FIRST EXECUTABLE STATEMENT
      IFAIL=0
      X = Y                                                             
      ISW = 1                                                           
      IF (X.GE.0.0D0) GO TO 5                                           
      ISW = -1                                                          
      X = -X                                                            
    5 IF (X.LT..477D0) GO TO 10                                         
      IF (X.LE.4.0D0) GO TO 30                                          
      IF (ISW .GT. 0) GO TO 40                                          
      IF (X.LT.XLARGE) GO TO 45                                         
      RES = 2.0D0                                                       
      GO TO 70                                                          
C                                  ABS(Y) .LT. .477, EVALUATE           
C                                  APPROXIMATION FOR ERFC               
   10 IF (X.LT.XMIN) GO TO 20                                           
      XSQ = X*X                                                         
      XNUM = P(4)*XSQ+P(5)                                              
      XDEN = XSQ+Q(4)                                                   
      DO 15 I = 1,3                                                     
         XNUM = XNUM*XSQ+P(I)                                           
         XDEN = XDEN*XSQ+Q(I)    
   15 CONTINUE                                                          
      RES = X*XNUM/XDEN                                                 
      GO TO 25                                                          
   20 RES = X*P(3)/Q(3)                                                 
   25 IF (ISW.EQ.-1) RES = -RES                                         
      RES = 1.0D0-RES                                                   
      GO TO 70                                                          
C                                  .477 .LE. ABS(Y) .LE. 4.0            
C                                  EVALUATE APPROXIMATION FOR ERFC      
   30 XSQ = X*X                                                         
      XNUM = P1(8)*X+P1(9)                                              
      XDEN = X+Q1(8)                                                    
      DO 35 I=1,7                                                       
         XNUM = XNUM*X+P1(I)                                            
         XDEN = XDEN*X+Q1(I)                                            
   35 CONTINUE                                                          
      RES = XNUM/XDEN                                                   
      GO TO 60                                                          
C                                  4.0 .LT. ABS(Y), EVALUATE            
C                                  MINIMAX APPROXIMATION FOR ERFC       
   40 IF (X.GT.XBIG) GO TO 65                                           
   45 XSQ = X*X                                                         
      XI = 1.0D0/XSQ                                                    
      XNUM= P2(5)*XI+P2(6)                                              
      XDEN = XI+Q2(5)                                                   
      DO 50 I = 1,4                                                     
         XNUM = XNUM*XI+P2(I)                                           
         XDEN = XDEN*XI+Q2(I)                                           
   50 CONTINUE                                                          
      RES = (SQRPI+XI*XNUM/XDEN)/X                                      
   60 RES = RES*DEXP(-XSQ)                                              
      IF (ISW.EQ.-1) RES = 2.0D0-RES                                    
      GO TO 70                                                          
   65 RES = 0.0D0                                                       
   70 S15ADF = RES                                                       
C
      RETURN                                                            
      END                                                               
