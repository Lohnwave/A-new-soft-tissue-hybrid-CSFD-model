SUBROUTINE Input(x,numd,nx,numnode,ndivx,ndivy,ndivz,ndivxq,ndivyq,ndivzq,&
nconn2,nquado,pAlf,Dmat,ALFs,numcell,numq,noCell,ncn,xc,&
npEBCnum,npEBC,pEBC,npNBCnum,npNBC,pNBC)
!------------------------------------------------------------------
! Input data from outside
! Output-all variables are output
!------------------------------------------------------------------
implicit real*8 (a-h,o-z)
common/para/xlength,ylength,zlength,p,young,anu,aimo
COMMON/rpim/ALFC,DC,Q,nRBF
common /basis/mbasis
!CHARACTER*100 NAM
dimension npEBC(4,1000),pEBC(3,1000),npNBC(4,1000),pNBC(3,1000)
dimension x(nx,numd),Dmat(6,6),noCell(8,ncn),xc(nx,numd)
!read(4,10)nam
read(4,*) xlength,ylength,zlength,young,anu,p
!read(4,10)nam
read(4,*)numnode,nconn2
!read(4,10)nam
read(4,*) ndivx,ndivy,ndivz
!read(4,10)nam
read(4,*)numq,numcell
!read(4,10)nam
read(4,*)ndivxq,ndivyq,ndivzq
!read(4,10)nam
read(4,*)nquado,pAlf
!read(4,10)nam
read(4,*)ALFs
!numgauss=nquado*nquado*nquado
!read(4,10)nam
do i=1,numnode
read(4,*)j,x(1,i),x(2,i),x(3,i)
enddo
!read(4,10)nam
do i=1,numq
read(4,*)j,xc(1,i),xc(2,i),xc(3,i)
enddo
!read(4,10)nam
do j=1,numcell
read(4,*)i,noCell(1,j),noCell(2,j),noCell(3,j),noCell(4,j),noCell(5,j),noCell(6,j),noCell(7,j),noCell(8,j)
enddo
!read(4,10)nam
read(4,*)npEBCnum
!read(4,10)nam
do i=1,npEBCnum
read(4,*)npEBC(1,i),npEBC(2,i),npEBC(3,i),npEBC(4,i),pEBC(1,i),pEBC(2,i),pEBC(3,i)
enddo
!read(4,10)nam
read(4,*)npNBCnum
!read(4,10)nam
do i=1,npNBCnum
read(4,*)npNBC(1,i),npNBC(2,i),npNBC(3,i),npNBC(4,i),pNBC(1,i),pNBC(2,i),pNBC(3,i)
enddo
!read(4,10)nam
READ(4,*)nRBF, alfc,dc, q
!read(4,10)nam
READ(4,*)mbasis
! ************* Compute material matrix D[] for the plane stress
you=young/(1.-anu*anu)
aimo=(1./12.)*xlength**3
Dmat(1,1)=you
Dmat(1,2)=anu*you
Dmat(1,3)=anu*you
Dmat(1,4)=0
Dmat(1,5)=0
Dmat(1,6)=0
Dmat(2,1)=anu*you
Dmat(2,2)=you
Dmat(2,3)=anu*you
Dmat(2,4)=0
Dmat(2,5)=0
Dmat(2,6)=0
Dmat(3,1)=anu*you
Dmat(3,2)=anu*you
Dmat(3,3)=you
Dmat(3,4)=0
Dmat(3,5)=0
Dmat(3,6)=0
Dmat(4,1)=0
Dmat(4,2)=0
Dmat(4,3)=0
Dmat(4,4)=0.5*(1.-anu)*you
Dmat(4,5)=0
Dmat(4,6)=0
Dmat(5,1)=0
Dmat(5,2)=0
Dmat(5,3)=0
Dmat(5,4)=0
Dmat(5,5)=0.5*(1.-anu)*you
Dmat(5,6)=0
Dmat(6,1)=0
Dmat(6,2)=0
Dmat(6,3)=0
Dmat(6,4)=0
Dmat(6,5)=0
Dmat(6,6)=0.5*(1.-anu)*you
!10 format(a40)
RETURN
END
