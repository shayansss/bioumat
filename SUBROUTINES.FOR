C
C     The following Fortran code was written by Seyed Shayan Sajjadinia. 
C     This code is the Fortran 77 version of the UMAT, FLOW, and SDVINI subroutines
C     of the cartilage model, I firstly proposed in my Master's thesis.
C     The model with minor modification was used in several publications.
C
C     If you use this code in your work, please cite the following article
C     as the reference, where you can also find more information about this code:
C
C     Sajjadinia, S. S., Haghpanahi, M., & Razi, M. (2019). Computational simulation of    
C     the multiphasic degeneration of the bone-cartilage unit during osteoarthritis via    
C     indentation and unconfined compression tests. Proceedings of the Institution of      
C     Mechanical Engineers, Part H: Journal of Engineering in Medicine.                    
C     https://doi.org/10.1177/0954411919854011                                             
C   
C
      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)

      DOUBLE PRECISION DEPTH,RPHI,G(10),RTHETA
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,TEN=10.D0,FOUR=4.D0,
     1 CONS1=5.235987755983D0,PI=3.14159265359D0,CONS2=2.61799387799D0)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     THIS SUBROUTINE CODE IS GENERATED TO CONTROL THE NON-HOMOGENEOUS AND VARIABLE MATERIAL PARAMETERS VIA STATEV
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C     MOST OF THE MATERIAL PARAMETERS ARE DEPENDING ON THE CARTLAGE DEGRADATION THAT ARE CONTROLED BY STATEV(1):
C     STATEV(1)=1 ==>> TOTALLY HEALTHY CARTILAGE
C     STATEV(1)=2 ==>> DEGRADED CARTILAGE WITHOUT FIBRILAR GEOMETRILAL ABNORMALITIES
C     STATEV(1)=3 ==>> DEGRADED CARTILAGE WITH FIBRILARATION IN DEEP ZONE
C     STATEV(1)=4 ==>> DEGRADED CARTILAGE WITH LOW FIBRILAR DISTORTION
C     STATEV(1)=5 ==>> DEGRADED CARTILAGE WITH HIGH FIBRILAR DISTORTION
C
      STATEV(1)=1
      DEPTH=-10D0 ! DEPTH OF THE CARTILAGE
      DEPTH=DBLE(COORDS(2)/DEPTH) ! NOMALIZATION OF DEPTH
C
C     VARIBLE PARAMETER DUE TO ANISOTROPIC AND NON-HOMOGENEOUS NATURE OF THE SOLID MATRIX:
C
      STATEV(2)=DBLE(1.4*(DEPTH**TWO)-1.1*DEPTH+0.59) ! DEPTH-DEPENDENT FIBER CONSTANT
      IF (STATEV(1).EQ.4) THEN
       CALL RANDOM-NUMBER(RPHI)
       RPHI=RPHI*PI/TWO
      ELSE
       IF (DEPTH.GT.0.3) THEN
        STATEV(4)=ONE
        STATEV(5)=ZERO
       ELSEIF (DEPTH.LE.0.3) THEN
        RPHI=DBLE(CONS1*DEPTH)
        STATEV(4)=DBLE(SIN(RPHI))
        IF (STATEV(1).EQ.5) THEN
         RTHETA=DBLE((PI/FOUR)-(CONS2*DEPTH))
        ELSEIF (STATEV(1).EQ.6) THEN
         RTHETA=DBLE((PI)-(FOUR*CONS2*DEPTH))
        ENDIF
        STATEV(5)=DBLE(COS(RPHI))
        STATEV(6)=ZERO
        IF (STATEV(1).EQ.(5.OR.6)) THEN
         STATEV(5)=DBLE(SIN(RTHETA)*STATEV(5))
         STATEV(6)=DBLE(COS(RTHETA)*STATEV(5))
        ENDIF
       ENDIF
       IF (STATEV(1).EQ.1) THEN ! DEPTH-DEPENDENT SOLID MATERIAL CONSTANT
        STATEV(3)=DBLE(0.1+0.2*DEPTH)
       ELSE
        STATEV(3)=DBLE(0.05+0.2*DEPTH)
       ENDIF
      ENDIF
      G(1)=0.005D0
      G(2)=0.01D0
      G(3)=0.025D0
      G(4)=0.035D0
      G(5)=0.042D0
      G(6)=0.048D0
      G(7)=0.053D0
      G(8)=0.058D0
      G(9)=0.06D0
      G(10)=0.06D0
      STATEV(7)=INT(DEPTH*9)+1
      STATEV(8)=G(STATEV(10)) ! GAG DEPTH-DEPENDENT MATERIAL CONSTANT
      STATEV(9)=DEPTH
C     THESE INITIALIZED STATE VARIBLES ARE ALSO USED FOR VERIFICATIONS
      DO i = 10, NSTATV
       STATEV(i)=ZERO
      ENDDO
C
C      TAFFETANI NEO-HOOKEAN MODEL (FOR VALIDATON)
C      G(1)=0.005D0
C      G(2)=0.008D0
C      G(3)=0.015D0
C      G(4)=0.05D0
C      G(5)=0.08D0
C      G(6)=0.1D0
C      G(7)=0.3D0
C      G(8)=0.65D0
C      G(9)=0.8D0
C      G(10)=0.8D0
C      STATEV(9)=G(INT(DEPTH*9)+1)
C
      RETURN
      END
C
C
C
C
      SUBROUTINE FLOW(H,SINK,U,KSTEP,KINC,TIME,NOEL,NPT,COORDS,
     1 JLTYP,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION TIME(2), COORDS(3)
      CHARACTER*80 SNAME
      H=1
      SINK=0
      IF ((COORDS(1).LE.1.25).AND.(COORDS(1).GE.-1.25)) THEN
       H=0
      ENDIF
      RETURN
      END
C
C
C
C
C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     LOCAL PARAMETERS AND VARIABLES
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      DOUBLE PRECISION TRANF(3,3),DET,IDENT(3,3),INDEX(2,6),ALPHA1,
     1 ALPHA2,IGAGD,NVEC0(27),RH,W1,W2,W3,W4,W5,W6,STRG(NTENS,NTENS),
     2 EPS,FV1(3),HH,NEWV1(3),NSTR(NTENS),STR,LANDA,DELTAV(NTENS),
     3 BVEC(NTENS),C,NS0,E1MP,E2MP,K1MP,CSTR,DFGRD(3,3),EP,STRS(NTENS),
     4 STATE(NSTATV),GAG,VV(NTENS)

      INTEGER i, j, k, l,r, m, n, K1, K2, K3, K4, K5, K6,KKK,DDS,FF
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1 TT=7.D0,SIX=6.D0,HALF=5.D-1,FFD=0.57735026919D0,
     2 R2P2=0.7071067811865D0,R3P3=0.5773502691896D0,TEN=10D0,
     3 CONST=0.4342944819D0)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     INITIALIZATION
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      DDS=1 ! THIS VARIABLE CONTROLS THE DDSDDE FORMULATION.
      ALPHA1=STATEV(8) ! DEPTH-DEPENDENT GAG CONSTANT
      ALPHA2=3.22D0 ! GAG CONSTANT
      C=3.009D0 ! FIBRILLAR RELATIVE DENSITY CONSTANT
      RH=STATEV(2) ! TOTAL FIBRILLAR DENSITY
      NS0=STATEV(3) ! DEPTH-DEPENDENT SOLID VOLUME FRACTION CONSTANT
      DO i = 1, 3
       DELTAV(i)=ONE ! KRONKER-DELTA IN VOIGT-NOTATION
       DO j = 1, 3
        IDENT(j,i) = ZERO ! KRONKER-DELTA 2ND ORDER TENSOR
        TRANF(j,i) = ZERO ! TRANSOPSE OF DEFORMATION GRADIENT TENSOR (DFGRD1) IN THE END OF THE INCREMNT
       ENDDO
       IDENT(i,i) = ONE
      ENDDO
      DO i=4,NTENS ! NTENS CONTROLS THE DIMENTIONALITY OF THE CODE
       DELTAV(i)=ZERO
      ENDDO
      DO i = 1,NTENS
       STRESS(i)=ZERO ! CAUCHY STRESS TENSOR THAT SHOULD BE UPDATED
       STRS(i)=ZERO ! FIBRILLAR STRESS TENSOR
       DO j = 1,NTENS
        DDSDDE(j,i)=ZERO ! JACOBIAN MATRIX TENSOR THAT SHOULD BE UPDATED.
       ENDDO
      ENDDO
      INDEX(1,1)=1 ! INDEX ARRAY ARE USED TO CIRCOMVENT ASSIGNING EQUAL COMPONENTS DUE TO SYMMYTRY OF HIYER ORDER TENSORS
      INDEX(2,1)=1
      INDEX(1,2)=2
      INDEX(2,2)=2
      INDEX(1,3)=3
      INDEX(2,3)=3
      INDEX(1,4)=1
      INDEX(2,4)=2
      INDEX(1,5)=1
      INDEX(2,5)=3
      INDEX(1,6)=2
      INDEX(2,6)=3
C
      NVEC0(1)=STATEV(5) ! NVEC ARRAY ARE THE DEFINED INITAIAL UNIT VECTOR DIRECTIONS ALONG FIBRILS.
      NVEC0(2)=STATEV(4)
      NVEC0(3)=STATEV(6)
C
      NVEC0(4)=-STATEV(5)
      NVEC0(5)=STATEV(4)
      NVEC0(6)=-STATEV(6)
C              
      NVEC0(7)=ONE
      NVEC0(8)=ZERO
      NVEC0(9)=ZERO
C              
      NVEC0(10)=ZERO
      NVEC0(11)=ONE
      NVEC0(12)=ZERO
C              
      NVEC0(13)=ZERO
      NVEC0(14)=ZERO
      NVEC0(15)=ONE
C              
      NVEC0(16)=FFD ! FFD SIGNIFIES 45 DEGREE
      NVEC0(17)=FFD
      NVEC0(18)=FFD
C              
      NVEC0(19)=-FFD
      NVEC0(20)=FFD
      NVEC0(21)=FFD
C              
      NVEC0(22)=FFD
      NVEC0(23)=-FFD
      NVEC0(24)=FFD
C              
      NVEC0(25)=FFD
      NVEC0(26)=FFD
      NVEC0(27)=-FFD
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     DDSDDE DERIVATION VIA PERTURBATION METHOD (INITIALIZATION)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      IF (DDS.EQ.1) THEN
       KKK=1
       EP=10D-8
       DO i = 1,3
        DO j = 1,3
         DFGRD(j,i)=DFGRD1(j,i)
        ENDDO
       ENDDO
       DO K6=1,NSTATV
        STATE(K6)=STATEV(K6)
       ENDDO
       i=1
       j=1
80     CONTINUE
       DO K1=1,3
        DO K2=1,3
         DFGRD1(K2,K1)=DFGRD(K2,K1)+(IDENT(K2,i)*DFGRD(j,K1)+IDENT(K2,j)
     1    *DFGRD(i,K1))*EP/TWO
        ENDDO
       ENDDO
90     CONTINUE
       DO K6=1,NTENS
        STRESS(K6)=ZERO
       ENDDO
       DO K6=1,NSTATV
        STATEV(K6)=STATE(K6)
       ENDDO
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     STRESS CALCULATIONS
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL TRANSPOSE(DFGRD1,TRANF)
      CALL VMATMUL(DFGRD1,TRANF,NTENS,BVEC) ! BVEC IS THE LEFT CAUCHY-GREEN OR FINGER DEFORMATION TENSOR
      CALL DETERMINANT(DFGRD1,DET) ! DET IS THE DETERMINANT OF DEFORMATION GRADIENT IN THE END OF THE INCREMENT.
C
      HH=NS0*(RH*C)/(TWO*C+TT) ! CNTRIBUTION OF OTHER CONSTANTS ON FIBRILAR STRESS
      E1MP=4.63D0 ! LINEAR MATERIAL CONSTANT OF FIBRILLAR PART
      E2MP=3670D0 ! NONLINEAR MATERIAL CONSTANT OF FIBRILLAR PART
      IF (STATEV(1).NE.1) THEN
       E2MP=3670/FOUR   ! FOR OA
      ENDIF
      E1MP=E1MP*HH
      E2MP=E2MP*HH
       DO i=0,8
       IF (i.EQ.2) THEN ! SECONDARY FIBRIL HAVE LOWER DENSITY BY C CONSTANT.
        E1MP=E1MP/C
        E2MP=E2MP/C
       ENDIF
       DO m = 1,3
        FV1(m)=ZERO ! FV VETOR IS THE INNER PRODUCT OF DFGRD1 AND NVEC0
        DO n = 1,3
         FV1(m)=DFGRD1(m,n)*NVEC0(3*i+n)+FV1(m)
        ENDDO
       ENDDO
       LANDA=SQRT(FV1(1)**TWO+FV1(2)**TWO+FV1(3)**TWO) ! LANDA IS THE ELONGATION
       EPS=LOG(LANDA) ! EPS IS THE FIBRIL LOGARITMIC STRIN
       IF (EPS.GT.ZERO) THEN
        DO n=1,3
         NEWV1(n)=FV1(n)/LANDA ! NEWV1 IS THE CURRENT FIBRIL DIRECTION
        ENDDO
        STR=(E1MP+E2MP*EPS)*EPS*LANDA/DET ! STR IS THE LOCAL FIBRIL STRESS
        DO K6=1,NTENS
         K3=INDEX(1,K6)
         K4=INDEX(2,K6)
         VV(K6)=NEWV1(K3)*NEWV1(K4) ! VV IS THE DYADIC PRODUCT OF CURRENT DIRECTION VECTORS THAT IS THE STRUCTRAL VECTOR
         STRS(K6)=STR*VV(K6) ! STRS IS THE GLOBAL FIBRIL STRESS
         STRESS(K6)=STRESS(K6)+STRS(K6)
        ENDDO
        IF (DDS.NE.1) THEN ! EXACT DDSDDE IMPLEMETATION FOR FIBRILLAR PART
         W1=(LANDA/DET)*E2MP*EPS
         W3=(ONE/EPS)-ONE
         DO K6=1,NTENS
          DO K5=1,NTENS
           DDSDDE(K5,K6)=(W1*VV(K5)+W3*STRS(K5))*VV(K6)+DDSDDE(K5,K6)
          ENDDO
         ENDDO
        ENDIF
       ENDIF
      ENDDO
      STATEV(10)=STRESS(2) ! S22 STRESS OF FIBRILLAR PART
C
      GM=0.723D0 ! GM IS THE NEO-HOOKEAN CONSTANT FOR PG CONRIBUTION
      GM=GM*NS0*(ONE-RH) ! CONTRIBUTION OF OTHER CONSTANTS
      W5=GM/DET
      W6=((LOG(DET)/SIX)*(((THREE*NS0/(DET-NS0))
     1 *((DET*LOG(DET)/(DET-NS0))-TWO))-FOUR)+(DET**(TWO/THREE)))*W5
      DO K6=1,NTENS
       NSTR(K6)=-DELTAV(K6)*W6+BVEC(K6)*W5
C      NSTR(K6)=STATEV(9)*(BVEC(K6)-DELTAV(K6))/DET ! TEFFANI NEO-HOOKEAN MODEL (FOR VALIDATION)
      ENDDO
      STATEV(11)=NSTR(2) ! S22 STRESS OF NON-FIBRILLAR PART
      GAG=ALPHA1*(DET**(-ALPHA2))
      STATEV(12)=GAG ! S22 STRESS OF GAG PART
      DO K6=1,3
       STRESS(K6)=NSTR(K6)-GAG+STRESS(K6)
      ENDDO
      DO K6=4,NTENS
       STRESS(K6)=NSTR(K6)+STRESS(K6)
      ENDDO
C
      IF (DDS.NE.1) THEN ! OTHER PARTS OF EXACT DDSDDE
       W3=(GM/TWO)*((FOUR*(DET**(TWO/THREE))-FOUR+(THREE*NS0/(DET-NS0))*
     1  (((DET*LOG(DET))/(DET-NS0))-TWO))/(THREE*DET)+((LOG(DET)-ONE)*
     2  NS0*LOG(DET))/((DET-NS0)**TWO))
       CALL TENF(STRESS,NTENS,STRG) ! STRG IS THE TENSOR FORM OF THE STRESS VECTOR
       W2=GAG*(ALPHA2-ONE)
       DO K6=1,NTENS
        K3 = INDEX(1,K6)
        K4 = INDEX(2,K6)
        DO K5=1,NTENS
         K1 = INDEX(1,K5)
         K2 = INDEX(2,K5)
         DDSDDE(K5,K6)=HALF*(IDENT(K4,K1)*STRG(K3,K2)+IDENT(K3,K2)*
     1    STRG(K1,K4)+IDENT(K3,K1)*STRG(K4,K2)+IDENT(K4,K2)*STRG(K1,K3))
     2    +W6*(IDENT(K1,K3)*IDENT(K2,K4)+IDENT(K1,K4)*IDENT(K2,K3))
     3    -W3*IDENT(K1,K2)*IDENT(K3,K4)+IDENT(K1,K2)*IDENT(K3,K4)*W2
     4    +(IDENT(K1,K4)*IDENT(K2,K3)+IDENT(K1,K3)*IDENT(K2,K4))*GAG
     5   +DDSDDE(K5,K6)
        ENDDO
       ENDDO
      ENDIF
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     DDSDDE DERIVATION VIA PERTURBATION METHOD (THE SECOND PART)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      IF (DDS.EQ.1) THEN
       IF (KKK.LT.NTENS) THEN
        DO K5=1,NTENS
         DDSDDE(K5,KKK)=STRESS(K5)*DET
        ENDDO
        KKK=KKK+1
        i=INDEX(1,KKK)
        j=INDEX(2,KKK)
        GO TO 80
       ENDIF
       IF (KKK.EQ.NTENS) THEN
        DO K5=1,NTENS
         DDSDDE(K5,KKK)=STRESS(K5)*DET
        ENDDO
        DO K1=1,3
         DO K2=1,3
          DFGRD1(K2,K1)=DFGRD(K2,K1)
         ENDDO
        ENDDO
        KKK=KKK+1
        GO TO 90
       ENDIF
       W2=ONE/EP
       W1=W2/DET
       DO K6=1,NTENS
        DO K5=1,NTENS
         DDSDDE(K5,K6)=W1*DDSDDE(K5,K6)-W2*STRESS(K5)
        ENDDO
       ENDDO
      ENDIF
C
      RETURN
      END
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     THE FOLLOWING SUBROUTINES ARE INTERNAL TO THE UMAT SUBROUTINE
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C     TRANSPOSE(A)
C
      SUBROUTINE TRANSPOSE(A,AT)
      INTEGER i , j
      DOUBLE PRECISION A(3,3), AT(3,3)
      Do i = 1 , 3
       DO j = 1 , 3
        AT(j,i) = A(i,j)
       ENDDO
      ENDDO
      RETURN
      END
C
C     DETERMINANT(A)
C
      SUBROUTINE DETERMINANT(A,DET)
      DOUBLE PRECISION A(3,3), DET
      DET=A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)
     1 -A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2)
     2 +A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2)
      RETURN
      END
C
C     INNER PRODUCT OF TWO MATRICES IN VOIGT-NOTATION
C
      SUBROUTINE VMATMUL(A,B,N,C)
      DOUBLE PRECISION A(3,3), B(3,3), C(N)
      INTEGER i , j
      C(1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      C(2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
      C(3)=A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
      C(4)=A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
      IF (N.EQ.6) THEN
      C(5)=A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
      C(6)=A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
      ENDIF
      RETURN
      END
C
C     MATRIX FORM OF A VECTOR
C
      SUBROUTINE TENF(V,N,M)
      DOUBLE PRECISION V(N),M(3,3)
      M(1,1)=V(1)
      M(2,1)=V(4)
      IF (N.EQ.6) THEN
      M(3,1)=V(5)
      ELSE
      M(3,1)=0
      ENDIF
      M(1,2)=V(4)
      M(2,2)=V(2)
      IF (N.EQ.6) THEN
      M(3,2)=V(6)
      M(1,3)=V(5)
      M(2,3)=V(6)
      ELSE
      M(3,2)=0
      M(1,3)=0
      M(2,3)=0
      ENDIF
      M(3,3)=V(3)
      RETURN
      END
C
      END
