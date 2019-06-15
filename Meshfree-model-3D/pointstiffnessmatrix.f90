SUBROUTINE PointStiffnessMatrix(ndex,weight,ajac,ph,Dmat,GSPk)
!----------------------------------------------------------------------------
! This subroutine to calculate sparse stiff matrix
! input--ndex: the number of nodes in the support domain;
! weight: weight of Gauss quadrature;
! ajac: Jacobian;
! dphix: first dirivetive of x of shape function;
! dphiy: first dirivetive of y of shape function;
! Dmat(3,3): the matrix of strain-stress;
! output--GSPk(2ndex,2ndex): sub-stiffness matrix of the Gauss point
!---------------------------------------------------------------------------
implicit real*8 (a-h,o-z)
dimension ph(10,ndex),Dmat(6,6),GSPk(3*ndex,3*ndex)
dimension bmat(6,3*ndex),dphix(ndex),dphiy(ndex),dphiz(ndex)
nb=3*ndex
do i=1,ndex
dphix(i)=ph(2,i)
dphiy(i)=ph(3,i)
dphiz(i)=ph(4,i)
enddo
do ib=1,6
do jb=1,nb
Bmat(ib,jb)=0.
enddo
enddo
do in=1,ndex
j=3*in-2
k=3*in-1
t=3*in
Bmat(1,j)=dphix(in)
Bmat(1,k)=0.
Bmat(1,t)=0.
Bmat(2,j)=0.
Bmat(2,k)=dphiy(in)
Bmat(2,t)=0.
Bmat(3,j)=0.
Bmat(3,k)=0.
Bmat(3,t)=dphiz(in)
Bmat(4,j)=0.
Bmat(4,k)=dphiz(in)
Bmat(4,t)=dphiy(in)
Bmat(5,j)=dphiz(in)
Bmat(5,k)=0.
Bmat(5,t)=dphix(in)
Bmat(6,j)=dphiy(in)
Bmat(6,k)=dphix(in)
Bmat(6,t)=0.
enddo
do ii=1,nb
do jj=1,nb
GSPk(ii,jj)=0.
enddo
enddo
do ii=1,nb
do jj=1,nb
do kk=1,6
do mm=1,6
GSPk(ii,jj)=GSPk(ii,jj)+weight*ajac*Bmat(kk,ii)* &
Dmat(kk,mm)*Bmat(mm,jj)
enddo
enddo
enddo
enddo
RETURN
END
