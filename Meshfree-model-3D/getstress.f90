SUBROUTINE GetStress(x,noCell,ds,Dmat,u2,alfs,nx,numnode,numgauss,&
xc,gauss,nquado,ng,numq,numcell, ENORM,Stressnode)
!----------------------------------------------------------------------------
! This subroutine to get the stress and energy error;
! input--numnode: total number of field nodes;
! numcell: number of cells;
! numq: total number of points for cells;
! alfs: coefficent of support support;
! x(nx,numnode): coordinates of all field nodes;
! xc(nx,numcell): coordinates of all points for cells;
! u2(2*numnode): displacement parameters;
! ds(nx,numnode): sizes of influence domain;
! Dmat(3,3): material matrix;
! nquado: number of Gauss points in a cell;
! gauss(nx,nquado): coefficients of Gauss points;
! numgauss: total number of Gauss points in all cells;
! output-- Enorm: energy error;
! Stressnode:stress for field nodes;
! compute out--Stress: stress for Gauss points;
! stressex, strne: exact stresse for beam problem.
!---------------------------------------------------------------------------
implicit real*8 (a-h,o-z)
common/para/xlength,ylength,p,young,anu,aimo
COMMON/rpim/ALFC,DC,Q,nRBF
common/basis/mbasis
dimension noCell(4,numcell),ds(nx,numnode),x(nx,numnode),u(2*numnode)
dimension xc(nx,numnode),gauss(nx,nquado)
dimension Dmat(3,3),u2(nx,numnode)
dimension Stressnode(3,numnode),strne(3,numnode)
dimension stress(3),stressex(3),err(3),Dinv(3,3),der(3)
integer, allocatable :: nv(:)
integer, allocatable :: ne(:)
real(8), allocatable :: ph(:,:)
real(8), allocatable :: gpos(:)
real(8), allocatable :: gs(:,:)
real(8), allocatable :: bmat(:,:)
allocate ( nv(1:numnode) )
allocate ( ne(1:2*numnode) )
allocate ( ph(1:10,1:3*numnode) )
allocate ( gpos(1:nx) )
allocate ( gs(1:ng,1:numgauss) )
allocate ( bmat(1:3,1:2*numnode) )
close(37)
open(37, file='midstr.dat')
do id=1,3
do jd=1,3
Dinv(id,jd)=Dmat(id,jd)
enddo
enddo
invd=3
ep=1.d-10
call GetINVASY(invd,invd,Dinv,EP)
do iu=1,numnode
ju=2*iu-1
ku=2*iu
u(ju)=u2(1,iu)
u(ku)=u2(2,iu)
enddo
enorm=0.
!****************Compute energy error
do 100 ibk=1,numcell
ind=0
call CellGaussPoints(ibk,numcell,nquado,numq,numgauss,&
xc,noCell,gauss,gs)
do 200 is=1,numgauss
do i=1,3
stress(i)=0.
stressex(i)=0.
enddo
ind=ind+1
gpos(1)= gs(1,is)
gpos(2)=gs(2,is)
weight=gs(3,is)
ajac=gs(4,is)
call SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)
do kph=1,3*ndex
do ik=1,10
ph(ik,kph)=0.
enddo
enddo
call RPIM_ShapeFunc_2D(gpos,x,nv,ph,nx,numnode,ndex,&
alfc,dc,q,nRBF, mbasis)
nb=2*ndex
do in=1,nb
ne(in)=0
enddo
do ine=1,ndex
n1=2*ine-1
n2=2*ine
ne(n1)=2*nv(ine)-1
ne(n2)=2*nv(ine)
enddo
do ib=1,3
do jb=1,nb
Bmat(ib,jb)=0.
enddo
enddo
do inn=1,ndex
j=2*inn-1
k=2*inn
Bmat(1,j)=ph(2,inn)
Bmat(1,k)=0.
Bmat(2,j)=0.
Bmat(2,k)=ph(3,inn)
Bmat(3,j)=ph(3,inn)
Bmat(3,k)=ph(2,inn)
enddo
do ii=1,3
do kk=1,3
do mm=1,nb
mn=ne(mm)
stress(ii)=stress(ii)+&
Dmat(ii,kk)*Bmat(kk,mm)*u(mn)
enddo
enddo
enddo
!****************Exact stress for beam problem
stressex(1)=(1./aimo)*p*(xlength-gpos(1))*gpos(2)
stressex(2)=0.
stressex(3)=-0.5*(p/aimo)*(0.25*ylength*ylength-gpos(2)*gpos(2))
do ier=1,3
err(ier)=stress(ier)-stressex(ier)
enddo
do jer=1,3
der(jer)=0.
do ker=1,3
der(jer)=der(jer)+Dinv(jer,ker)*err(ker)
enddo
enddo
err2=0.
do mer=1,3
err2=err2+weight*ajac*(0.5*der(mer)*err(mer))
enddo
enorm=enorm+err2
200 continue
100 continue
!****************Compute nodal stresses
write(2,*)'Stress of field nodes'
do is=1,numnode
gpos(1)= x(1,is)
gpos(2)=x(2,is)
do ii=1,3
Stressnode(ii,is)=0.
enddo
ndex=0
call SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)
do kph=1,3*ndex
do ik=1,10
ph(ik,kph)=0.
enddo
enddo
call RPIM_ShapeFunc_2D(gpos,x,nv,ph,nx,numnode,ndex,&
alfc,dc,q,nRBF, mbasis)
nb=2*ndex
do in=1,nb
ne(in)=0
enddo
do ine=1,ndex
n1=2*ine-1
n2=2*ine
ne(n1)=2*nv(ine)-1
ne(n2)=2*nv(ine)
enddo
do ib=1,3
do jb=1,nb
Bmat(ib,jb)=0.
enddo
enddo
do inn=1,ndex
j=2*inn-1
k=2*inn
Bmat(1,j)=ph(2,inn)
Bmat(1,k)=0.
Bmat(2,j)=0.
Bmat(2,k)=ph(3,inn)
Bmat(3,j)=ph(3,inn)
Bmat(3,k)=ph(2,inn)
enddo
do ii=1,3
do kk=1,3
do mm=1,nb
mn=ne(mm)
Stressnode(ii,is)=Stressnode(ii,is)+&
Dmat(ii,kk)*Bmat(kk,mm)*u(mn)
enddo
enddo
enddo
strne(1,is)=(1./aimo)*p*(xlength-gpos(1))*gpos(2)
strne(2,is)=0.
strne(3,is)=-0.5*(p/aimo)*(0.25*ylength*ylength-gpos(2)*gpos(2))
write(2,220) is,Stressnode(1,is),Stressnode(2,is),Stressnode(3,is)
if(abs(gpos(1)-24).le.1.d-5) then
write(37,240) is,gpos(2),Stressnode(1,is),Stressnode(2,is), &
Stressnode(3,is),strne(1,is),strne(2,is),strne(3,is)
endif
enddo
enorm=dsqrt(enorm)
write(2,230) enorm
220 format(1x,i5,3e20.5)
230 format(1x,'Energy Error=',e20.10)
240 format(1x,i5,f8.3,6e15.5)
deallocate ( nv)
deallocate ( ne)
deallocate ( ph)
deallocate ( gpos)
deallocate ( gs )
deallocate ( bmat)
RETURN
END
