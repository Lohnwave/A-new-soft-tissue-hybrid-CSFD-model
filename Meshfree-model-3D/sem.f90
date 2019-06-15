
subroutine sme(n,a,b,x)


!--------------------------------------program comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  高斯列主元消去法
!    
!-----------------------------------------------------
!  In put data  files :
!       1.  fin.txt  输入方程系数
!       2.
!  Output data files  :
!       1. fout.txt  计算结果
!       2.
!-----------------------------------------------------
!  Post Script :
!       1.    需要准备输入数据
!     
!-----------------------------------------------------
use m_gauss

implicit real*8(a-h,o-z)



integer::i,j
real*8::A(N,N),b(N),x(N)

call solve(A,b,x,N)
return
end 