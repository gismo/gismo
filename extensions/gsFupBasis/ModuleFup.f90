MODULE FUP_0_16_D
!   
!  U PROGRAMU KOJI UKLJUCUJE (MODULE FUP_0_16) POTREBNOJE DA PRVA IZVRSNA NAREDBA BUDE
! ******************  CALL  RACUN  ******************
!   
!
!  MODUL DONOSI U GLAVNI PROGRAM POLJA S VRIJEDNOSTIMA DERIVACIJA FUNKCIJE
!  FUPn(x),  NA RAZMAKU 2**(-M), M = 16,  n = 0-0,...,16-16
!         !!!  N A J B R Z A   V A R I J A N T A  !!!
!  Derivacije za pojedinu funkciju Fupn(x) se izracunavaju od nultog do n-tog reda




   PUBLIC RACUN, FUPN  !

   PRIVATE UPTURBOX, FUP00, FUP01, FUP02, FUP03, FUP04, FUP05, FUP06, FUP07,  &
                      FUP08, FUP09, FUP10, FUP11, FUP12, FUP13, FUP14, FUP15, FUP16, &
                      NFUP, VERTEX, XPOINT, DELTAX, KOD, K,M,N,KK,J00,J01,J02,J03,  &
                      WORK, D,C00,DX,DX0 


   CONTAINS

   SUBROUTINE RACUN
    
!
!  PODPROGRAM VRACA POLJA S VRIJEDNOSTIMA DERIVACIJA FUNKCIJE FUPn(x) NA RAZMAKU 2**(-M), M = 16
!  n = 0-0,...,16-16         !!!  N A J B R Z A   V A R I J A N T A  !!!
!
   INTEGER(4) :: K,M,N,KK,J00,J01,J02,J03

   real(kind=8)  D,C00,DX,DX0
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16

   DIMENSION D(32)
   real(8), ALLOCATABLE :: WORK(:,:)
   !(-2*65536:65536,0:16)

!   PARAMETER ( M = 16 )
!   PARAMETER ( DX = 1.0D0 )  ! DEFAULT, INACE SE MORA ZADATI KONKRETNA VRIJEDNOST ' DX '
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)

!   CALL CPU_TIME(T0)
!   WRITE(*,*) T0

   M  = 16
   DX = 1.D0   ! DEFAULT, INACE SE MORA ZADATI KONKRETNA VRIJEDNOST ' DX '


!!! Calculate coefficient Delta (+-1)

   D(1) = 1.D0
   DO K = 1,16
   D(2*K-1) = D(K)
   D(2*K) = -D(K)
   END DO

!!! Fill array of U (value of UP function - 2E16)

          CALL UPTURBOX(M)

!IZRACUNAVANJE SVIH DERIVACIJA FUNKCIJE up(x)


      ALLOCATE (WORK(-2*65536:65536,0:16))
         K = -2**M
      DO K = -2**M,2**M
      WORK(K,0) = FUP_00(K, 0)
     
      END DO
         N = 1
      DO N = 1,M  !!!+1  !!!!
              K = 0
           DO K = 0,2**(M-N+1) 
           WORK(-2**M+K,N) = 2.0D0**((N*(N+1))/2)*FUP_00(-2**M+K*2**N,0)
           END DO
      END DO
         N = 2
      DO N = 2,M   !!!+1   !!!
           DO K = 0,2**(M-N+1)
                           KK = 2
                        DO KK = 2,N
     WORK(-2**M+2**(M-N+1)*(KK-1)+K, N) = D(KK)*WORK(-2**M+K,N)
                         END DO
           END DO
      END DO


!    Izracunavanje vrijednosti nulte i prvih (N+1) derivacija funkcije Fupn(x)
           N = 0
        DO N = 0,M
              J00 = 2**M
              J01 = ((N+2)*2**(M-N))/2
  J02 = 0
  J03 = 2**(M-N)
 C00 = 2.0D0**((N*(N+1))/2) 

  DX0 = 2.0D0**(-N)
                      
                           IF(N == 1) THEN
                    KK = 0
                 DO KK = 0,N+1  !!!
                     K = 0
                 DO  K = 0,J01            
 FUP_01(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-WORK(-J00-J03+K,KK))
     FUP_01( J01-K,KK) = (-1.0D0)**KK*FUP_01(-J01+K,KK)
                 END DO
                 END DO
          
                      ELSE IF(N == 2) THEN

                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
  DO  K = 0,J01           
FUP_02(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-2.0D0*WORK(-J00-J03+K,KK))
FUP_02( J01-K,KK) = (-1.0D0)**KK*FUP_02(-J01+K,KK)
                 END DO
                 END DO

                      ELSE IF(N == 3) THEN
 			
                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
  DO  K = 0,J01          
FUP_03(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-3.0D0*WORK(-J00-J03+K,KK)+ &
              4.0D0*WORK(-J00-2*J03+K,KK))
 FUP_03( J01-K,KK) = (-1.0D0)**KK*FUP_03(-J01+K,KK)
                 END DO
                 END DO
 
                      ELSE IF(N == 4) THEN

                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
     DO  K = 0,J01           
 FUP_04(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-4.0D0*WORK(-J00-J03+K,KK)+ &
                  7.0D0*WORK(-J00-2*J03+K,KK))
     FUP_04( J01-K,KK) = (-1.0D0)**KK*FUP_04(-J01+K,KK)
                 END DO
                 END DO

                      ELSE IF(N == 5) THEN

                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
     DO  K = 0,J01   
 FUP_05(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-5.0D0*WORK(-J00-J03+K,KK)+ &
                   11.0D0*WORK(-J00-2*J03+K,KK)-15.0D0*WORK(-J00-3*J03+K,KK))
     FUP_05( J01-K,KK) = (-1.0D0)**KK*FUP_05(-J01+K,KK)
                 END DO
                 END DO

                      ELSE IF(N == 6) THEN

                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
     DO  K = 0,J01              
 FUP_06(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-6.0D0*WORK(-J00-J03+K,KK)+ &
                 16.0D0*WORK(-J00-2*J03+K,KK)-26.0D0*WORK(-J00-3*J03+K,KK))
     FUP_06( J01-K,KK) = (-1.0D0)**KK*FUP_06(-J01+K,KK)
                 END DO
                 END DO
 
                     ELSE IF(N == 7) THEN

                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
    DO  K = 0,J01              
 FUP_07(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-7.0D0*WORK(-J00-J03+K,KK)+ &
  22.0D0*WORK(-J00-2*J03+K,KK)-42.0D0*WORK(-J00-3*J03+K,KK)+ &
 58.0D0*WORK(-J00-4*J03+K,KK))
     FUP_07( J01-K,KK) = (-1.0D0)**KK*FUP_07(-J01+K,KK)
                 END DO
                 END DO

                      ELSE IF(N == 8) THEN

                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
  DO  K = 0,J01             
FUP_08(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-8.0D0*WORK(-J00-J03+K,KK)+ &
 29.0D0*WORK(-J00-2*J03+K,KK)-64.0D0*WORK(-J00-3*J03+K,KK)+ &
 100.0D0*WORK(-J00-4*J03+K,KK))
     FUP_08( J01-K,KK) = (-1.0D0)**KK*FUP_08(-J01+K,KK)
                 END DO
                 END DO

                      ELSE IF(N == 9) THEN

                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
 DO  K = 0,J01              
FUP_09(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-9.0D0*WORK(-J00-J03+K,KK)+ &
 37.0D0*WORK(-J00-2*J03+K,KK)-93.0D0*WORK(-J00-3*J03+K,KK)+ &
164.0D0*WORK(-J00-4*J03+K,KK)-228.0D0*WORK(-J00-5*J03+K,KK))
     FUP_09( J01-K,KK) = (-1.0D0)**KK*FUP_09(-J01+K,KK)
                 END DO
                 END DO

                      ELSE IF(N == 10) THEN
                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
 DO  K = 0,J01           
FUP_10(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-10.0D0*WORK(-J00-J03+K,KK)+ &
           46.0D0*WORK(-J00-2*J03+K,KK)-130.0D0*WORK(-J00-3*J03+K,KK)+ &
  257.0D0*WORK(-J00-4*J03+K,KK)-392.0D0*WORK(-J00-5*J03+K,KK))
     FUP_10( J01-K,KK) = (-1.0D0)**KK*FUP_10(-J01+K,KK)
                 END DO
                 END DO

 ELSE IF(N == 11) THEN

                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
   DO  K = 0,J01               
FUP_11(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-11.0D0*WORK(-J00-J03+K,KK)+ &
 56.0D0*WORK(-J00-2*J03+K,KK)-176.0D0*WORK(-J00-3*J03+K,KK)+ &
  387.0D0*WORK(-J00-4*J03+K,KK)-649.0D0*WORK(-J00-5*J03+K,KK)+ &
 904.0D0*WORK(-J00-6*J03+K,KK))
     FUP_11( J01-K,KK) = (-1.0D0)**KK*FUP_11(-J01+K,KK)
                 END DO
                 END DO

                      ELSE IF(N == 12) THEN

                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
  DO  K = 0,J01              
 FUP_12(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-12.0D0*WORK(-J00-J03+K,KK)+ &
  67.0D0*WORK(-J00-2*J03+K,KK)-232.0D0*WORK(-J00-3*J03+K,KK)+ &
 563.0D0*WORK(-J00-4*J03+K,KK)-1036.0D0*WORK(-J00-5*J03+K,KK)+ &
   1553.0D0*WORK(-J00-6*J03+K,KK))
     FUP_12( J01-K,KK) = (-1.0D0)**KK*FUP_12(-J01+K,KK)

                 END DO
                 END DO

                      ELSE IF(N == 13) THEN

                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
  DO  K = 0,J01              
 FUP_13(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-13.0D0*WORK(-J00-J03+K,KK)+ &
                  79.0D0*WORK(-J00-2*J03+K,KK)-299.0D0*WORK(-J00-3*J03+K,KK)+ &
  795.0D0*WORK(-J00-4*J03+K,KK)-1599.0D0*WORK(-J00-5*J03+K,KK)+ &
   2589.0D0*WORK(-J00-6*J03+K,KK)-3601.0D0*WORK(-J00-7*J03+K,KK))
     FUP_13( J01-K,KK) = (-1.0D0)**KK*FUP_13(-J01+K,KK)
                 END DO
                 END DO

                      ELSE IF(N == 14) THEN

                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
  DO  K = 0,J01              
FUP_14(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-14.0D0*WORK(-J00-J03+K,KK)+ &
92.0D0*WORK(-J00-2*J03+K,KK)-378.0D0*WORK(-J00-3*J03+K,KK)+ &
1094.0D0*WORK(-J00-4*J03+K,KK)-2394.0D0*WORK(-J00-5*J03+K,KK)+ &
4188.0D0*WORK(-J00-6*J03+K,KK)-6190.0D0*WORK(-J00-7*J03+K,KK))
     FUP_14( J01-K,KK) = (-1.0D0)**KK*FUP_14(-J01+K,KK)

                 END DO
                 END DO

                      ELSE IF(N == 15) THEN

                    KK = 0
                 DO KK = 0,N+1 !!!
                     K = 0
 DO  K = 0,J01               
FUP_15(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-15.0D0*WORK(-J00-J03+K,KK)+ &
106.0D0*WORK(-J00-2*J03+K,KK)-470.0D0*WORK(-J00-3*J03+K,KK)+ &
1472.0D0*WORK(-J00-4*J03+K,KK)-3488.0D0*WORK(-J00-5*J03+K,KK)+ &
6582.0D0*WORK(-J00-6*J03+K,KK)-10378.0D0*WORK(-J00-7*J03+K,KK)+ &
14384.0D0*WORK(-J00-8*J03+K,KK))
     FUP_15( J01-K,KK) = (-1.0D0)**KK*FUP_15(-J01+K,KK)
                 END DO
                 END DO

                      ELSE IF(N == 16) THEN

                    KK = 0
                 DO KK = 0,N !!!
                     K = 0
DO  K = 0,J01             
FUP_16(-J01+K,KK) = (DX0**KK)*C00*(WORK(-J00+K,KK)-16.0D0*WORK(-J00-J03+K,KK)+ &
121.0D0*WORK(-J00-2*J03+K,KK)-576.0D0*WORK(-J00-3*J03+K,KK)+ &
1942.0D0*WORK(-J00-4*J03+K,KK)-4960.0D0*WORK(-J00-5*J03+K,KK)+ &
10070.0D0*WORK(-J00-6*J03+K,KK)-16960.0D0*WORK(-J00-7*J03+K,KK)+ &
24762.0D0*WORK(-J00-8*J03+K,KK))
     FUP_16( J01-K,KK) = (-1.0D0)**KK*FUP_16(-J01+K,KK)
                 END DO
                 END DO

                      END IF
    
      END DO
!
      DEALLOCATE (WORK)

!
      END SUBROUTINE RACUN
   
!!!  Calculate value of UP function in the arbitrary points

      SUBROUTINE UPTURBOX(M)
  
       
real(8) UN,FAK,UNN,UN0,SUMAK,FUP_00,ZERO
COMMON FUP_00(0:131072, 0:1)
      DIMENSION UN(0:20),FAK(0:20),UNN(0:20),UN0(0:20)
      INTEGER(4) M,N,L,I,K
      
DATA UN0/1.0D0, 1.0D0, 5.0D0, 1.0D0, 143.0D0, 19.0D0, 1153.0D0,&
      583.0D0,1616353.0D0,132809.0D0, 134926369.0D0, 46840699.0D0,&
      67545496213157.0D0,4068990560161.0D0,411124285571171.0D0,&
      1204567303451311.0D0,73419800947733963069.0D0,&
      4146897304424408411.0D0,86773346866163284480799923.0D0,&
      18814360006695807527868793.0D0,&
      539741515875650532056045666422369.0D0/
!!!
      DATA UNN/1.0D0, 2.0D0, 72.0D0, 288.0D0, 2073600.0D0,&
      33177600.0D0, 561842749440.0D0, 179789679820800.0D0,&
      704200217922109440000.0D0, 180275255788060016640000.0D0,&
      1246394851358539387238350848000.0D0,&
      6381541638955721662660356341760000.0D0,&
      292214732887898713986916575925267070976000000.0D0,&
      1196911545908833132490410294989893922717696000000.0D0,&
      17524030168305511965050671651660013242599473361715200000.0D0,&
     15791254065263462941946461238743871133171237435708801024000000.0D0,&
626048168100066478643636623385103067560649311175417997219699097600000000.0D0,&
48488455061807039788823800613832680933046479303954410932297509162188800000000.0D0,  &
2924907327984663493179931480281060829152039746976389598046631610524565901849531514880000000.0D0,  &
3833734532936058133780799789833992049986161537156893373951680984546759018872217947183513600000000.0D0,  &
1391026453346497029228605426710671587340398341822365591105864134203768975246595308345805415518935449600000000000.0D0/
!!!
      DATA FAK/1.0D0,1.0D0,2.0D0,6.0D0,24.0D0,120.0D0,720.0D0,5040.0D0, &
          40320.0D0,362880.0D0,3628800.0D0,39916800.0D0,479001600.0D0,  &
      6227020800.0D0,87178291200.0D0,1307674368000.0D0,&
      20922789888000.0D0,355687428096000.0D0,6402373705728000.0D0,&
      121645100408832000.0D0,2432902008176640000.0D0/
!!!
      DATA ZERO/0.0D0/

           I = 0
        DO I = 0,M
                    UN(I) = UN0(I)/UNN(I)
FAK(I)= 2.0D0**((I*(I+1))/2)/FAK(I)
        END DO 

                  FUP_00  = 0.0D0
 FUP_00(1,0) = UN(M)
FUP_00(2,0) = UN(M-1)

 N = 1
 DO N = 1,M 
K = 2**N                      
 DO K = 2**N,2**(N+1)
                          SUMAK = ZERO
                                         L = 0
                                      DO L = 0,M-N 
           SUMAK = SUMAK + FAK(L)*UN(M-N-L)*(2.0D0**(-M)*(K-2**N))**L
                                      END DO
                                      FUP_00(K, 0) = SUMAK - FUP_00(K-2**N, 0)
                          END DO
                  END DO
		  
		  
		  
DO K = 0, 65536
                         FUP_00(K,1) = 2.0*FUP_00(2*K,0)
                         FUP_00(K+65536,1) = -FUP_00(K,1)
END DO 
 
     END SUBROUTINE UPTURBOX



    real(8) FUNCTION FUPN(NFUP, VERTEX, XPOINT, DELTAX, KOD)
	 
                      
    INTEGER(4) ::  NFUP,KOD
real(8)    ::  VERTEX, XPOINT, DELTAX

    SELECT CASE (NFUP)

    CASE ( 0)
    FUPN = FUP00(VERTEX, XPOINT, DELTAX, KOD)
    CASE ( 1)
    FUPN = FUP01(VERTEX, XPOINT, DELTAX, KOD)
CASE ( 2)
FUPN = FUP02(VERTEX, XPOINT, DELTAX, KOD)
    CASE ( 3)    
FUPN = FUP03(VERTEX, XPOINT, DELTAX, KOD)
    CASE ( 4)    
FUPN = FUP04(VERTEX, XPOINT, DELTAX, KOD)
    CASE ( 5)    
FUPN = FUP05(VERTEX, XPOINT, DELTAX, KOD)
    CASE ( 6)    
FUPN = FUP06(VERTEX, XPOINT, DELTAX, KOD)
    CASE ( 7)    
FUPN = FUP07(VERTEX, XPOINT, DELTAX, KOD)
    CASE ( 8)    
FUPN = FUP08(VERTEX, XPOINT, DELTAX, KOD)
    CASE ( 9)    
FUPN = FUP09(VERTEX, XPOINT, DELTAX, KOD)
    CASE (10)    
FUPN = FUP10(VERTEX, XPOINT, DELTAX, KOD)
    CASE (11)    
FUPN = FUP11(VERTEX, XPOINT, DELTAX, KOD)
    CASE (12)    
FUPN = FUP12(VERTEX, XPOINT, DELTAX, KOD)
    CASE (13)    
FUPN = FUP13(VERTEX, XPOINT, DELTAX, KOD)
    CASE (14)    
FUPN = FUP14(VERTEX, XPOINT, DELTAX, KOD)
    CASE (15)    
FUPN = FUP15(VERTEX, XPOINT, DELTAX, KOD)
    CASE (16)    
FUPN = FUP16(VERTEX, XPOINT, DELTAX, KOD)
    END SELECT

    END FUNCTION FUPN



    real(8) FUNCTION FUP00(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, INDEX, DVANF, KOD
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD=0, 1  -  indeks reda trazene derivacije 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 0
   DVANF = 65536     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP00 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP00 = FUP_00(INDEX, KOD)
   END IF
   END FUNCTION FUP00




    real(8) FUNCTION FUP01(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 1
   DVANF = 32768     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP01 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP01 = FUP_01(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP01




    real(8) FUNCTION FUP02(VERTEX, XPOINT, DELTAX, KOD)
	 
INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 2
   DVANF = 16384     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY) .GE.  float((DVANF*(NFUP+2)/2))) THEN
              FUP02 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP02 = FUP_02(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP02




    real(8) FUNCTION FUP03(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2,3
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 3
   DVANF = 8192     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP03 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP03 = FUP_03(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP03





    real(8) FUNCTION FUP04(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2,3,4
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 4
   DVANF = 4096     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP04 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP04 = FUP_04(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP04


    real(8) FUNCTION FUP05(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2,3,4,5
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 5
   DVANF = 2048     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP05 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP05 = FUP_05(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP05


    real(8) FUNCTION FUP06(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2,3,4,5,6
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 6
   DVANF = 1024     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP06 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP06 = FUP_06(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP06


    real(8) FUNCTION FUP07(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2,3,4
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 7
   DVANF = 512     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP07 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP07 = FUP_07(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP07


    real(8) FUNCTION FUP08(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2, ... ,7,8
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 8
   DVANF = 256     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP08 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP08 = FUP_08(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP08


    real(8) FUNCTION FUP09(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2, ... ,9
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 9
   DVANF = 128     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP09 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP09 = FUP_09(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP09


    real(8) FUNCTION FUP10(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2, ... ,10
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 10
   DVANF = 64     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP10 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP10 = FUP_10(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP10


    real(8) FUNCTION FUP11(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2, ... ,11
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 11
   DVANF = 32     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP11 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP11 = FUP_11(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP11

    real(8) FUNCTION FUP12(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1, ... ,12
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 12
   DVANF = 16     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP12 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP12 = FUP_12(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP12


    real(8) FUNCTION FUP13(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1, ... ,13
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 13
   DVANF =  8     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP13 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP13 = FUP_13(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP13


    real(8) FUNCTION FUP14(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2, ... ,14
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 14
   DVANF =  4     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP14 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP14 = FUP_14(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP14


    real(8) FUNCTION FUP15(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2, ... ,15
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 15
   DVANF =  2   !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP15 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP15 = FUP_15(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP15


    real(8) FUNCTION FUP16(VERTEX, XPOINT, DELTAX, KOD)
	 
    INTEGER(4) NFUP, KOD, INDEX, DVANF
real(8) VERTEX, XPOINT, DELTAX,  DUMMY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Funkcijski podprogram za izracunavanje vrijednosti funkcije i prvih n derivacija        
!   za bazne funkcije Fupn(x) reda n = 0,1,2, ... ,16
!     
!           NFUP   -  red "n" odabrane bazne funkcije Fupn(x)
!           VERTEX -  koordinata tjemena odabrane bazne funkcije Fupn(x)
!           XPOINT -  koordinata tocke u kojoj se izracunava vrijednost derivacije reda <=n
!           DELTAX -  duljina karakteristicnog odsjecka
!           KOD    -  indeks reda trazene derivacije 0,1,2, ... ,16
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(8)  FUP_00, FUP_01, FUP_02, FUP_03, FUP_04, FUP_05, FUP_06, FUP_07, FUP_08, &
            FUP_09, FUP_10, FUP_11, FUP_12, FUP_13, FUP_14, FUP_15, FUP_16
   COMMON   FUP_00(-65536:65536,0: 1), &
            FUP_01(-49152:49152,0: 2), FUP_02(-32768:32768,0: 3), FUP_03(-20480:20480,0: 4), &
            FUP_04(-12288:12288,0: 5), FUP_05( -7168: 7168,0: 6), FUP_06( -4096: 4096,0: 7), &
            FUP_07( -2304: 2304,0: 8), FUP_08( -1280: 1280,0: 9), FUP_09(  -704:  704,0:10), &
            FUP_10(  -384:  384,0:11), FUP_11(  -208:  208,0:12), FUP_12(  -112:  112,0:13), &
            FUP_13(   -60:   60,0:14), FUP_14(   -32:   32,0:15), FUP_15(   -17:   17,0:16), &
            FUP_16(    -9:    9,0:17)
   NFUP  = 16
   DVANF =  1     !  2**(16-NFUP)
   DUMMY = (XPOINT-VERTEX)/DELTAX* float(DVANF)

   IF(abs(DUMMY).GE. float((DVANF*(NFUP+2)/2))) THEN
              FUP16 = 0.0D0
   ELSE
       INDEX = Nint(DUMMY)
       FUP16 = FUP_16(INDEX,KOD)/DELTAX**KOD
   END IF
   END FUNCTION FUP16



END MODULE FUP_0_16_D    
