
subroutine sme(n,a,b,x)


!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  ��˹����Ԫ��ȥ��
!    
!-----------------------------------------------------
!  In put data  files :
!       1.  fin.txt  ���뷽��ϵ��
!       2.
!  Output data files  :
!       1. fout.txt  ������
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.    ��Ҫ׼����������
!     
!-----------------------------------------------------
use m_gauss

implicit real*8(a-h,o-z)



integer::i,j
real*8::A(N,N),b(N),x(N)

call solve(A,b,x,N)
return
end 