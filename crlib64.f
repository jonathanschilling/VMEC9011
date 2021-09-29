      subroutine spaxpy(n, sa, sx, sy, ind)

      real sx(*), sy(*), sa
      integer n, ind(*)
      integer i

      if(n.le.0) return
      do i = 1, n
        sy(ind(i)) = sa*sx(i) + sy(ind(i))
      enddo
      end
      SUBROUTINE SOLR(NP,A,INCA,B,INCB,C,INCC)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION A(1),B(1),C(1)
      IA = NP * IABS(INCA)
      IF (INCA .GT. 0) IA = 1
      IB = NP * IABS(INCB)
      IF (INCB .GT. 0) IB = 1
      IC = NP * IABS(INCC)
      IF (INCC .GT. 0) IC = 1
      IC1 = IC + INCC
      IC2 = IC1 + INCC
      DO 10 K = 1,NP-2
         C(IC2) = A(IA) * C(IC1) + B(IB) * C(IC)
         IA = IA + INCA
         IB = IB + INCB
         IC = IC + INCC
         IC1 = IC1 + INCC
         IC2 = IC2 + INCC
   10 CONTINUE
      RETURN
      END
      SUBROUTINE MXV (A, NA, B, NB, C)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION A(NA,NB), B(NB   ), C(NA   )
         DO 20 I = 1, NA
            C( I ) = 0.0
            DO 10 K = 1, NB
   10       C( I ) = C( I ) + A(I,K) * B( K )
   20    CONTINUE
      RETURN
      END
      SUBROUTINE MXM (A, NA, B, NB, C, NC)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION A(NA,NB), B(NB,NC), C(NA,NC)
      DO 30 J = 1, NC
         DO 20 I = 1, NA
            C(I,J) = 0.0
            DO 10 K = 1, NB
   10       C(I,J) = C(I,J) + A(I,K) * B(K,J)
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
      FUNCTION CVMGT (A, B, C)
      LOGICAL C
      CVMGT = B
      IF ( C ) CVMGT = A
      RETURN
      END
      FUNCTION CVMGM (A, B, C)
      CVMGM = B
      IF (C .LT. 0.0) CVMGM = A
      RETURN
      END
      FUNCTION CVMGP (A, B, C)
      CVMGP = B
      IF (C .GE. 0.0) CVMGP = A
      RETURN
      END
      FUNCTION CVMGZ (A, B, C)
      CVMGZ = B
      IF (C .EQ. 0.0) CVMGZ = A
      RETURN
      END
      FUNCTION CVMGN (A, B, C)
      CVMGN = B
      IF (C .NE. 0.0) CVMGN = A
      RETURN
      END
      SUBROUTINE ORDERS (MODE, S1, E, IX, NP, I1, I8, I2)
C
C     ALTERNATIVUM TO CRAY-SORT ROUTINE
C     SORTS FLOATING POINT NUMBERS ACCORDING TO THEIR SIZE
C
C     E(I) - KEYS TO BE SORTED
C     IX(I) - INDEX : E(IX(1)) ON RETURN IS THE SMALLEST
C                     E(IX(2)) ON RETURN IS THE NEXT TO THE SMALLEST
C
C     NP NUMBER OF KEYS TO BE SORTED
C     S1 IS A SCRATCH STORAGE OF LENGTH 256 FOR THE CRAY-VERSION
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION IX(1), E(1)
      IF (MODE .NE. 2) PRINT *,' MODE NE 2 ', MODE
      IF (I1   .NE. 1) PRINT *,' I1   NE 1 ', I1
      IF (I8   .NE. 8) PRINT *,' I8   NE 8 ', I8
      IF (I2   .NE. 1) PRINT *,' I2   NE 1 ', I2
      DO 10 I = 1, NP
   10 IX(I) = I
      DO 100 I = 1, NP
         DO 50 J = I+1, NP
            IF (E(IX(I)) .LT. E(IX(J))) GO TO 50
            IT = IX(I)
            IX(I) = IX(J)
            IX(J) = IT
   50    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE MINV (AB, N, ND, SCR, DET, EPS, M, MODE)
C
C     LINPACK SIMULATION OF CRAY-MINV
C
C     AB      AUGMENTED MATRIX  FIRST DIMENSION IS ND
C     N       ORDER OF MATRIX
C     ND      FIRST DIMENSION OF MATRIX AB
C     SCR     SCRATCH VECTOR OF DIMENSION AT LEAST 2*N
C     DET     DETERMINANT
C     EPS     TOLERANCE FOR THE PRODUCT OF PIVOT ELEMENTS
C                         ----> NOT NEEDED HERE
C     M       > 0 NUMBER OF SYSTEMS OF LINEAR EQUATIONS TO SOLVE
C             = 0 DETERMINANT IS COMPUTED
C     MODE    = 1 AB IS OVERWRITTEN WITH INV (AB)
C             = 0 DETERMINANT AS WELL AS INVERSE ARE NOT COMPUTED
C                 ONLY EQUATIONS ARE SOLVED  (IF ANY)
C                 MATRIX IS NOT SAVED IN THIS CASE
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION AB(ND,1), SCR(1), DT(2)
C
C     PERFORM L-U DECOMPOSITION
C
C     CALL SGECO (AB(1,1), ND, N, SCR, RCOND, SCR(N+1))
      CALL SGEFAL(AB(1,1), ND, N, SCR, INFO)
C
C     SOLVE LINEAR EQUATIONS IF ANY
C
      IF (M .EQ. 0) GO TO 20
      DO 10 I = 1,M
         CALL SGESLL(AB(1,1), ND, N, SCR, AB(1,N+I), 0)
   10 CONTINUE
   20 IF (MODE .EQ. 0) RETURN
C
C     COMPUTE DETERMINANT AND INVERSE
C
      CALL SGEDIL(AB, ND, N, SCR(1), DT, SCR(N+1), 11)
      DET = DT(1) * (10.0**DT(2))
      RETURN
      END
      SUBROUTINE SGEFAL(A, LDA, N, IPVT, INFO)
      IMPLICIT double precision (A-H,O-Z)
      INTEGER LDA,N,IPVT(1),INFO
      DIMENSION A(LDA,1)
C
C     SGEFA FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
C
C     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF RCOND IS NOT NEEDED
C     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA)
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE MATRIX TO BE FACTORED
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF ARRAY A
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX A
C
C     ON RETURN
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT
C                THE FACTORIZATION CAN BE WRITTEN A = L*U WHERE
C                L IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND U IS UPPER TRIANGULAR
C
C        IPVT    INTEGER(N)
C                AN INTEGER OF PIVOT INDICES
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE
C                = K  IF  U(K,K) .EQ. 0.0   THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE RCOND IN SGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY
C
C     LINPACK  THIS VERSION DATED 08/14/78
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY, SSCAL, ISAMAX
C
C     INTERNAL VARIABLES
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = ISAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0E0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0E0/A(K,K)
            CALL SSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL SAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
      SUBROUTINE SGESLL(A,LDA,N,IPVT,B,JOB)
      IMPLICIT double precision (A-H,O-Z)
      INTEGER LDA,N,IPVT(1),JOB
      DIMENSION A(LDA,1),B(1)
C
C     SGESL SOLVES THE REAL SYSTEM
C     A * X = B OR TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY SGECO OR SGEFA
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE OUPUT FROM SGECO OR SGEFA
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF ARRAY A
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX A
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SGECO OR SGEFA
C
C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR
C
C        JOB     INTEGER
C                = 0         TO SOLVE A*X = B
C                = NONZERO   TO SOLVE TRANS(A)*X = B WHERE
C                            TRANS(A) IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR X
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL. TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA. IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
C        OR SGEFA HAS SET INFO .EQ. 0
C
C     TO COMPUTE INVERSE(A) * C WHERE C IS A MATRIX
C     WITH P COLUMNS
C          CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
C          IF (RCOND IS TOO SMALL) GO TO ...
C          DO 10 J = 1, P
C             CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
C       10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SDOT
C
C     INTERNAL VARABLES
C
      INTEGER K,KB,L,NM1
C
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE A * X = B
C        FIRST SOLVE L * Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL SAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL SAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE TRANS(A) * X = B
C        FIRST SOLVE TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = SDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + SDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE SGEDIL (A,LDA,N,IPVT,DET,WORK,JOB)
      IMPLICIT double precision (A-H,O-Z)
      INTEGER LDA,N,IPVT(1),JOB
      DIMENSION A(LDA,1),DET(2),WORK(1)
C
C     SGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY SGECO OR SGEFA.
C
C     ON ENTRY
C
C        A       READ(LDA, N)
C                THE OUTPUT FROM SGECO OR SGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY A
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX A
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SGECO OR SGEFA
C
C        WORK    REAL(N)
C                WORK VECTOR. CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY
C                = 10   DETERMINANT ONLY
C     ON RETURN
C
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED
C
C        DET     REAL(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH 1.0 .LE. ABS(DET(1) .LT. 10.0
C                OR DET(1) .EQ. 0.0
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUROUTINES ARE CALLED CORRECTLY
C        AND IF SGECO HAS SET RCOND .GT. 0.0 OR SGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SSCAL,SSWAP
C     FORTRAN ABS,MOD
C
C     INTERNAL VARIABLES
C
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0E0
         DET(2) = 0.0E0
         TEN = 10.0E0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0E0) GO TO 60
   10       IF (ABS(DET(1)) .GE. 1.0E0) GO TO 20
               DET(1) = TEN * DET(1)
               DET(2) = DET(2) - 1.0E0
            GO TO 10
   20       CONTINUE
   30       IF (ABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0E0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0E0/A(K,K)
            T = -A(K,K)
            CALL SSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0E0
               CALL SAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0E0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL SAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL SSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
      SUBROUTINE SGBFAL (ABD,LDA,N,ML,MU,IPVT,INFO)
C
C     SGBFA FACTORS A REAL BAND MATRIX BY ELIMINATION.
C
C     SGBFA IS USUALLY CALLED BY SGBCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C
C     ON ENTRY
C
C        ABD     REAL(LDA, N)
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
C                SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C                0 .LE. ML .LT. N .
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. MU .LT. N .
C                MORE EFFICIENT IF  ML .LE. MU .
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF
C                     CALLED.  USE  RCOND  IN SGBCO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     BAND STORAGE
C
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
C           WILL SET UP THE INPUT.
C
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SSCAL,ISAMAX
C     FORTRAN MAX0,MIN0
C
      IMPLICIT double precision (A-H,O-Z)
      INTEGER LDA,N,ML,MU,IPVT(1),INFO
      DIMENSION ABD(LDA,1)
C
C
      M = ML + MU + 1
      INFO = 0
C
C     ZERO INITIAL FILL-IN COLUMNS
C
      J0 = MU + 2
      J1 = MIN0(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0E0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
C
C        ZERO NEXT FILL-IN COLUMN
C
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0E0
   40       CONTINUE
   50    CONTINUE
C
C        FIND L = PIVOT INDEX
C
         LM = MIN0(ML,N-K)
         L = ISAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (ABD(L,K) .EQ. 0.0E0) GO TO 100
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0E0/ABD(M,K)
            CALL SSCAL(LM,T,ABD(M+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL SAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0E0) INFO = N
      RETURN
      END
      SUBROUTINE SGBSLL(ABD,LDA,N,ML,MU,IPVT,B,JOB)
C
C     SGBSL SOLVES THE REAL BAND SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY SGBCO OR SGBFA.
C
C     ON ENTRY
C
C        ABD     REAL(LDA, N)
C                THE OUTPUT FROM SGBCO OR SGBFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SGBCO OR SGBFA.
C
C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SGBCO HAS SET RCOND .GT. 0.0
C        OR SGBFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SDOT
C     FORTRAN MIN0
C
      IMPLICIT double precision (A-H,O-Z)
      INTEGER LDA,N,ML,MU,IPVT(1),JOB
      DIMENSION ABD(LDA,1),B(1)
C
      INTEGER K,KB,L,LA,LB,LM,M,NM1
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE L*Y = B
C
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN0(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL SAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL SAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = SDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML,N-K)
               B(K) = B(K) + SDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE SPBFAL(ABD,LDA,N,M,INFO)
      IMPLICIT double precision (A-H,O-Z)
      INTEGER LDA,N,M,INFO
      DIMENSION ABD(LDA,*)
C
C     SPBFA FACTORS A REAL SYMMETRIC POSITIVE DEFINITE MATRIX
C     STORED IN BAND FORM.
C
C     SPBFA IS USUALLY CALLED BY SPBCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C
C     ON ENTRY
C
C        ABD     REAL(LDA, N)
C                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER
C                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE
C                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE
C                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. M + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        M       INTEGER
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. M .LT. N .
C
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND
C                FORM, SO THAT  A = TRANS(R)*R .
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT
C                     POSITIVE DEFINITE.
C
C     BAND STORAGE
C
C           IF  A  IS A SYMMETRIC POSITIVE DEFINITE BAND MATRIX,
C           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT.
C
C                   M = (BAND WIDTH ABOVE DIAGONAL)
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-M)
C                      DO 10 I = I1, J
C                         K = I-J+M+1
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 03/11/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LABS.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SDOT
C     FORTRAN MAX0,SQRT
C
C     INTERNAL VARIABLES
C
      INTEGER IK,J,JK,K,MU
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
         DO 30 J = 1, N
            INFO = J
            S = 0.0E0
            IK = M + 1
            JK = MAX0(J-M,1)
            MU = MAX0(M+2-J,1)
            DO 10 K = MU, M
               T = ABD(K,J) - SDOT(K-MU,ABD(IK,JK),1,ABD(MU,J),1)
               T = T/ABD(M+1,JK)
               ABD(K,J) = T
               S = S + T*T
               IK = IK - 1
               JK = JK + 1
   10       CONTINUE
            S = ABD(M+1,J) - S
C     ......EXIT
            IF (S .LE. 0.0E0) GO TO 40
            ABD(M+1,J) = SQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
      SUBROUTINE SPBSLL(ABD,LDA,N,M,B)
      IMPLICIT double precision (A-H,O-Z)
      INTEGER LDA,N,M
      DIMENSION ABD(LDA,*),B(*)
C
C     SPBSL SOLVES THE REAL SYMMETRIC POSITIVE DEFINITE BAND
C     SYSTEM  A*X = B
C     USING THE FACTORS COMPUTED BY SPBCO OR SPBFA.
C
C     ON ENTRY
C
C        ABD     REAL(LDA, N)
C                THE OUTPUT FROM SPBCO OR SPBFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        M       INTEGER
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C
C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
C        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
C        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
C        CORRECTLY AND  INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SPBCO(ABD,LDA,N,RCOND,Z,INFO)
C           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL SPBSL(ABD,LDA,N,C(1,J))
C        10 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 03/11/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LABS.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SDOT
C     FORTRAN MIN0
C
C     INTERNAL VARIABLES
C
      INTEGER K,LA,LB,LM
C
C     SOLVE TRANS(R)*Y = B
C
      DO 10 K = 1, N
         LM = MIN0(K-1,M)
         LA = M + 1 - LM
         LB = K - LM
         T = SDOT(LM,ABD(LA,K),1,B(LB),1)
         B(K) = (B(K) - T)/ABD(M+1,K)
   10 CONTINUE
C
C     SOLVE R*X = Y
C
      DO 20 K = N, 1, -1
         LM = MIN0(K-1,M)
         LA = M + 1 - LM
         LB = K - LM
         B(K) = B(K)/ABD(M+1,K)
         T = -B(K)
         CALL SAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   20 CONTINUE
      RETURN
      END
      FUNCTION SSUM(IM,A,INCA)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION A(1)
      SSUM = 0.0
      JA = IM * IABS(INCA)
      IF (INCA .GT. 0) JA = 1
      DO 10 J = 1,IM
         SSUM = SSUM + A(JA)
         JA = JA + INCA
   10 CONTINUE
      RETURN
      END
      SUBROUTINE SCOPY(IM,A,INCA,B,INCB)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION A(1),B(1)
      JA = IM * IABS(INCA)
      IF (INCA .GT. 0) JA = 1
      JB = IM * IABS(INCB)
      IF (INCB .GT. 0) JB = 1
      DO 10 J = 1,IM
         B(JB) = A(JA)
         JA = JA + INCA
         JB = JB + INCB
   10 CONTINUE
      RETURN
      END
      INTEGER FUNCTION ISAMAX(N,SX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK 3/11/78.
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION SX(1)
      INTEGER I, INCX,IX,N
C
      ISAMAX = 0
      IF (N .LT. 1) RETURN
      ISAMAX = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C     CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      SMAX = ABS(SX(1))
      IX = IX + INCX
      DO 10 I = 2, N
         IF (ABS(SX(IX)) .LE. SMAX) GO TO 5
         ISAMAX = I
         SMAX = ABS(SX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     CODE FOR INCREMENT EQUAL 1
C
   20 SMAX = ABS(SX(1))
      DO 30 I = 2, N
         IF (ABS(SX(I)) .LE. SMAX) GO TO 30
         ISAMAX = I
         SMAX = ABS(SX(I))
   30 CONTINUE
      RETURN
      END
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE
C     JACK DONGARRA, LINPACK, 3/11/78
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION SX(1),SY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C     PRINT *, ' SAXPY-N-A ', N, SA
C
      IF (N .LE. 0) RETURN
      IF (SA .EQ. 0.0) RETURN
      IF ((INCX .EQ. 1) .AND. (INCY .EQ. 1)) GO TO 20
C
C     CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C     NOT EQUAL TO ONE
C
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1, N
         SY(IY) = SY(IY) + SA * SX(IX)
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     CODE FOR BOTH INCREMENTS EQUAL TO ONE
C
C
C     CLEAN UP LOOP
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1, M
         SY(I) = SY(I) + SA * SX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1, N, 4
         SY(I  ) = SY(I  ) + SA * SX(I  )
         SY(I+1) = SY(I+1) + SA * SX(I+1)
         SY(I+2) = SY(I+2) + SA * SX(I+2)
         SY(I+3) = SY(I+3) + SA * SX(I+3)
   50 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION SDOT(N,SX,INCX,SY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION SX(*),SY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      STEMP = 0.0E0
      SDOT = 0.0E0
      IF (N .LE. 0) RETURN
      IF ((INCX .EQ. 1) .AND. (INCY .EQ. 1)) GO TO 20
C
C     CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C     NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1, N
         STEMP = STEMP + SX(IX) * SY(IY)
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      SDOT = STEMP
      RETURN
C
C     CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C     CLEAN UP LOOP
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1, M
         STEMP = STEMP + SX(I) * SY(I)
   30 CONTINUE
      IF (N .LT. 5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         STEMP = STEMP + SX(I)*SY(I) + SX(I+1)*SY(I+1) + SX(I+2)*SY(I+2)
     1                               + SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4)
   50 CONTINUE
   60 SDOT = STEMP
      RETURN
      END
      SUBROUTINE SSCAL (N,SA,SX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO 1
C     JACK DONGARRA, LINPACK, 3/11/78
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION SX(1)
      INTEGER I,INCX,M,MP1,N,NINCX
C
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C     CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1, NINCX, INCX
         SX(I) = SA*SX(I)
   10 CONTINUE
      RETURN
C
C     CODE FOR INCREMENT EQUAL TO 1
C
C
C     CLEAN UP LOOP
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1, M
         SX(I) = SA*SX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         SX(I  ) = SA*SX(I  )
         SX(I+1) = SA*SX(I+1)
         SX(I+2) = SA*SX(I+2)
         SX(I+3) = SA*SX(I+3)
         SX(I+4) = SA*SX(I+4)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE SSWAP (N,SX,INCX,SY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION SX(1), SY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF (N .LT. 0) RETURN
      IF ((INCX .EQ. 1) .AND. (INCY .EQ. 1)) GO TO 20
C
C     CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C     TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1, N
         STEMP = SX(IX)
         SX(IX) = SY(IY)
         SY(IY) = STEMP
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C     CLEAN UP LOOP
C
   20 M = MOD(N,3)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1, M
         STEMP = SX(I)
         SX(I) = SY(I)
         SY(I) = STEMP
   30 CONTINUE
      IF (N .LT. 3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
         STEMP   = SX(I  )
         SX(I  ) = SY(I  )
         SY(I  ) = STEMP
         STEMP   = SX(I+1)
         SX(I+1) = SY(I+1)
         SY(I+1) = STEMP
         STEMP   = SX(I+2)
         SX(I+2) = SY(I+2)
         SY(I+2) = STEMP
   50 CONTINUE
      RETURN
      END
      INTEGER FUNCTION ISMAX(N,SX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. VALUE.
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION SX(1)
      INTEGER I, INCX,IX,N
C
      ISMAX = 0
      IF (N .LT. 1) RETURN
      ISMAX = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C     CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      SMAX = SX(1)
      IX = IX + INCX
      DO 10 I = 2, N
         IF (SX(IX) .LE. SMAX) GO TO 5
         ISMAX = I
         SMAX = SX(IX)
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     CODE FOR INCREMENT EQUAL 1
C
   20 SMAX = SX(1)
      DO 30 I = 2, N
         IF (SX(I) .LE. SMAX) GO TO 30
         ISMAX = I
         SMAX = SX(I)
   30 CONTINUE
      RETURN
      END
      INTEGER FUNCTION ISMIN(N,SX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MIN. VALUE.
C
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION SX(1)
      INTEGER I, INCX,IX,N
C
      ISMIN = 0
      IF (N .LT. 1) RETURN
      ISMIN = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C     CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      SMIN = SX(1)
      IX = IX + INCX
      DO 10 I = 2, N
         IF (SX(IX) .GE. SMIN) GO TO 5
         ISMIN = I
         SMIN = SX(IX)
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     CODE FOR INCREMENT EQUAL 1
C
   20 SMIN = SX(1)
      DO 30 I = 2, N
         IF (SX(I) .GE. SMIN) GO TO 30
         ISMIN = I
         SMIN = SX(I)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE   RS   (XXX)
      IMPLICIT double precision (A-H,O-Z)
      XXX = 999999
      RETURN
      END
      SUBROUTINE SECOND (III)
      III = 999999
      RETURN
      END
      FUNCTION   RANF()
      IMPLICIT double precision (A-H,O-Z)
      RANF= 999999
      RETURN
      END
