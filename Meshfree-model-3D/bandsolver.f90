SUBROUTINE BANDSOLVER(A,F,B,N,L,IL,nmat)
!------------------------------------------------------------------
! Slover for banded linear equations
!------------------------------------------------------------------
implicit real*8 (a-h,o-z)
DIMENSION A(N,N),F(N)
DIMENSION B(N,nmat),d(n,1)
M=1
LP1=L+1
DO I=1,N
DO K=1,IL
B(I,K)=0.
IF(I.LE.LP1) B(I,K)=A(I,K)
IF(I.GT.LP1.AND.I.LE.(N-L)) B(I,K)=A(I,I+K-LP1)
IF(I.GT.(N-L).AND.(I+K-LP1).LE.N) B(I,K)=A(I,I+K-LP1)
ENDDO
ENDDO
DO I=1,N
D(I,1)=F(I)
ENDDO
IT=1
IF (IL.NE.2*L+1) THEN
IT=-1
WRITE(*,*)'***FAIL***'
RETURN
END IF
LS=L+1
DO 100 K=1,N-1
P=0.0
DO I=K,LS
IF (ABS(B(I,1)).GT.P) THEN
P=ABS(B(I,1))
IS=I
END IF
ENDDO
IF (P+1.0.EQ.1.0) THEN
IT=0
WRITE(*,*)'***FAIL***'
RETURN
END IF
DO J=1,M
T=D(K,J)
D(K,J)=D(IS,J)
D(IS,J)=T
ENDDO
DO J=1,IL
T=B(K,J)
B(K,J)=B(IS,J)
B(IS,J)=T
ENDDO
DO J=1,M
D(K,J)=D(K,J)/B(K,1)
ENDDO
DO J=2,IL
B(K,J)=B(K,J)/B(K,1)
ENDDO
DO I=K+1,LS
T=B(I,1)
DO J=1,M
D(I,J)=D(I,J)-T*D(K,J)
ENDDO
DO J=2,IL
B(I,J-1)=B(I,J)-T*B(K,J)
ENDDO
B(I,IL)=0.0
ENDDO
IF (LS.NE.N) LS=LS+1
100 CONTINUE
IF (ABS(B(N,1))+1.0.EQ.1.0) THEN
IT=0
WRITE(*,*)'***FAIL***'
RETURN
END IF
DO J=1,M
D(N,J)=D(N,J)/B(N,1)
ENDDO
JS=2
DO 150 I=N-1,1,-1
DO K=1,M
DO J=2,JS
D(I,K)=D(I,K)-B(I,J)*D(I+J-1,K)
ENDDO
ENDDO
IF (JS.NE.IL) JS=JS+1
150 CONTINUE
if(it.le.0) write(*,*) "BandSolver failed"
DO I=1,N
F(I)=D(I,1)
ENDDO
RETURN
END