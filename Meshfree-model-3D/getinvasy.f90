SUBROUTINE GetInvasy(N,MA,A,EPS)
!----------------------------------------------------------------------------
! This subroutine to get INVARSION OF A(N,N) USING THE GAUSS-JODON METHOD.
! MATRIX A MUST BE DEFINITE BUT MAY BE ASYMMETRIC.
! input--N: dimension of A;
! MA: max number of rows of A;
! EPS: tolerance;
! Input and output--A[N,N]: the matrix in the input and invarsion in output;
!---------------------------------------------------------------------------
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION A(MA,N)
DO 10 K=1,N
C=A(K,K)
IF(DABS(C).LE.EPS)pause
C=1.0/C
A(K,K)=1.0
DO J=1,N
A(K,J)=A(K,J)*C
ENDDO
DO 10 I=1,N
IF(I.EQ.K)GOTO 10
C=A(I,K)
A(I,K)=0.0
DO J=1,N
A(I,J)=A(I,J)-A(K,J)*C
ENDDO
10 CONTINUE
RETURN
END
