SUBROUTINE CellGaussPoints(ibk,numcell,k,numq,numgauss,xc,noCell,gauss,gs)
!----------------------------------------------------------------------------
! This subroutine to set up Gauss points,Jacobian and weights for a cell
! input--ibk: the No. of the consider cell;
! numq: number of points for background cells;
! numcell: number of background cells;
! numgauss: number of Gauss points in a cell;
! k: number of Gauss points used, numgauss=k*k for 2D cell;
! xc(nx,numq): coordinates of points for background cells;
! noCell(ng,numcell): No. of points to form this cell;
! gauss(2,k): coefficients of Gauss points;
! nx,ng: parameters are defined in file parameter.h.
! output--gs(ng,numgauss): coordinate of the Gauss points, weight and Jacobian
!---------------------------------------------------------------------------
implicit real*8 (a-h,o-z)
include 'parameter.h'
dimension xc(nx,numq),noCell(ng,numcell),gauss(2,k),gs(ng,numgauss)
dimension psiJ(ng),etaJ(ng),JetJ(ng),xe(ng),ye(ng),ze(ng),aN(ng),aNJpsi(ng),aNJeta(ng),aNJJet(ng)
index=0
psiJ(1)=-1.
psiJ(2)=1.
psiJ(3)=1.
psiJ(4)=-1.
psiJ(5)=-1.
psiJ(6)=1.
psiJ(7)=1.
psiJ(8)=-1.


etaJ(1)=-1.
etaJ(2)=-1.
etaJ(3)=1.
etaJ(4)=1.
etaJ(5)=-1.
etaJ(6)=-1.
etaJ(7)=1.
etaJ(8)=1.

JetJ(1)=-1.
JetJ(2)=-1.
JetJ(3)=-1.
JetJ(4)=-1.
JetJ(5)=1.
JetJ(6)=1.
JetJ(7)=1.
JetJ(8)=1.


l=k
ie=ibk
do j=1,ng
je=noCell(j,ie)
xe(j)=xc(1,je)
ye(j)=xc(2,je)
ze(j)=xc(3,je)
enddo
do 10 i=1,l
do 10 j=1,l
do 10 t=1,l
index=index+1
eta=gauss(1,i)
psi=gauss(1,j)
Jet=gauss(1,t)
do ik=1,ng
aN(ik)=1./8*(1.+psi*psiJ(ik))*(1.+eta*etaJ(ik))*(1+Jet*JetJ(ik))
aNJpsi(ik)=1./8*psiJ(ik)*(1.+eta*etaJ(ik))*(1+Jet*JetJ(ik))
aNJeta(ik)=1./8*etaJ(ik)*(1.+psi*psiJ(ik))*(1+Jet*JetJ(ik))
aNJJet(ik)=1./8*JetJ(ik)*(1.+psi*psiJ(ik))*(1.+eta*etaJ(ik))
enddo
xpsi=0.
ypsi=0.
zpsi=0.
xeta=0.
yeta=0.
zeta=0.
xJet=0.
yJet=0.
zJet=0.

do jk=1,ng
xpsi=xpsi+aNJpsi(jk)*xe(jk)
ypsi=ypsi+aNJpsi(jk)*ye(jk)
zpsi=zpsi+aNJpsi(jk)*ze(jk)
xeta=xeta+aNJeta(jk)*xe(jk)
yeta=yeta+aNJeta(jk)*ye(jk)
zeta=zeta+aNJeta(jk)*ze(jk)
xJet=xJet+aNJJet(jk)*xe(jk)
yJet=yJet+aNJJet(jk)*ye(jk)
zJet=zJet+aNJJet(jk)*ze(jk)
enddo
ajcob=xpsi*(yeta*zJet-zeta*yJet)-ypsi*(xeta*zJet-zeta*xJet)+zpsi*(xeta*yJet-yeta*xJet)
xq=0.
yq=0.
zq=0.
do kk=1,ng
xq=xq+aN(kk)*xe(kk)
yq=yq+aN(kk)*ye(kk)
zq=zq+aN(kk)*ze(kk)
enddo
gs(1,index)=xq
gs(2,index)=yq
gs(3,index)=zq
gs(4,index)=gauss(2,i)*gauss(2,j)*gauss(2,k)
gs(5,index)=ajcob
10 continue
RETURN
END
