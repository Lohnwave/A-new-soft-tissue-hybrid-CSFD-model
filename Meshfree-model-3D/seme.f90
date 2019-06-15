module m_gauss
!----------------------------------------module coment
!  Version     :  V1.0    
!  Coded by    :  syz 
!  Date        :  2010-4-8
!-----------------------------------------------------
!  Description : ��˹����Ԫ��ȥ��ģ��
!    
!-----------------------------------------------------
!  Contains    :
!      1.   solve  ��������
!      2.
!-----------------------------------------------------

contains

subroutine solve(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  ��˹����Ԫ��ȥ��
!                 Ax=b
!-----------------------------------------------------
!  Input  parameters  :
!       1.   A(N,N)ϵ������
!       2.   b(N)������
!       3.   N����ά��
!  Output parameters  :
!       1.  x  ���̵ĸ�
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)

integer::i,k,N
integer::id_max  !��Ԫ�ر��

real*8::A(N,N),b(N),x(N)

real*8::Aup(N,N),bup(N)

!AbΪ�������  [Ab]
real*8::Ab(N,N+1)

real*8::vtemp1(N+1),vtemp2(N+1)

Ab(1:N,1:N)=A

Ab(:,N+1)=b


!##########################################################
!  ����� ����Ԫ��ȥ���ĺ��Ĳ���
do k=1,N-1

    elmax=dabs(Ab(k,k))
    id_max=k
    
    !���Ϊ������Ԫ��	
    !��γ������ҪĿ�Ĳ���Ϊ�˸�ֵ���Ԫ�ظ�elmax������Ϊ���ҳ����Ԫ�ض�Ӧ�ı��

	
	do i=k+1,n
      if (dabs(Ab(i,k))>elmax) then
         elmax=Ab(i,k)

         id_max=i
      end if          
    end do

    
 !���ˣ��Ѿ���ɲ������Ԫ�أ���������Ժ���  ��k�н��� 
 !��������Ԫ�أ���������
    vtemp1=Ab(k,:)
    vtemp2=Ab(id_max,:)
   
    
    Ab(k,:)=vtemp2
    Ab(id_max,:)=vtemp1   
!
!����һ�����Ϊ��������Ԫ�أ���������Ժ󼴰�����Ԫ������
!#########################################################
  
   do i=k+1,N
  
     temp=Ab(i,k)/Ab(k,k)
     
     Ab(i,:)=Ab(i,:)-temp*Ab(k,:)
   
   end do

end do

!-----------------------------
! ������һ����Ab�Ѿ���Ϊ������ʽ�ľ���
!            | *  *  *  *  # |
!     [A b]= | 0  *  *  *  # |
!            | 0  0  *  *  # |
!            | 0  0  0  *  # |
!
Aup(:,:)=Ab(1:N,1:N)

bup(:)=Ab(:,N+1)

!�����������Ƿ�����Ļش�����
call uptri(Aup,bup,x,n)

end subroutine solve



subroutine uptri(A,b,x,N)
!---------------------------------subroutine  comment
!  Version   :  V1.0    
!  Coded by  :  syz 
!  Date      :  2010-4-8
!-----------------------------------------------------
!  Purpose   :  �����Ƿ�����Ļش�����
!                 Ax=b
!-----------------------------------------------------
!  Input  parameters  :
!       1.   A(N,N)ϵ������
!       2.   b(N)������
!       3.   N����ά��
!  Output parameters  :
!       1.  x  ���̵ĸ�
!       2.
!  Common parameters  :
!
!----------------------------------------------------

implicit real*8(a-z)

integer::i,j,N

real*8::A(N,N),b(N),x(N)

x(N)=b(N)/A(N,N)

!�ش�����
do i=n-1,1,-1
   
    x(i)=b(i)
   do j=i+1,N
    x(i)=x(i)-a(i,j)*x(j)
   end do
    x(i)=x(i)/A(i,i)

end do

end subroutine uptri

end module m_gauss
