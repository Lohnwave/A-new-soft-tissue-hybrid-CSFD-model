SUBROUTINE naturalBC_distributed(numnode,numq,in,jn,alfs,x,xc,ds, &
gauss,nquado,force,Fi)
!----------------------------------------------------------------------------
! This subroutine to enforce point natural bc's;
! input-numnode, numq, in,jn,alfs,x,xc,ds,gauss, nquado.
! input and output-- force{}:force vector.
!---------------------------------------------------------------------------
implicit real*8 (a-h,o-z)
include 'parameter.h'
common/para/xlength,ylength,p,young,anu,aimo
COMMON/rpim/ALFC,DC,Q,nRBF
COMMON /basis/mbasis
dimension xc(nx,numq),gauss(2,nquado),force(2*numnode),x(nx,numnode)
dimension ph(10,numnode),gpos(2),nv(numnode),ds(nx,numnode)
ax=0.5*(xc(1,in)-xc(1,jn))
ay=0.5*(xc(2,in)-xc(2,jn))
bx=0.5*(xc(1,in)+xc(1,jn))
by=0.5*(xc(2,in)+xc(2,jn))
do il=1,nquado
gpos(1)=ax*gauss(1,il)+bx
gpos(2)=ay*gauss(1,il)+by
weight=gauss(2,il)
ajac=0.5*sqrt((xc(1,in)-xc(1,jn))**2+(xc(2,in)-xc(2,jn))**2)
ximo=(1./12.)*xlength**3
ty=(-Fi/(2.*ximo))*(xlength*xlength/4.-(gpos(1)-24)*(gpos(1)-24))
call SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)
do kph=1,ndex
do ik=1,10
ph(ik,kph)=0.
enddo
enddo
call RPIM_ShapeFunc_2D(gpos,x,nv,ph,nx,numnode,ndex,&
alfc,dc,q,nRBF, mbasis)
do ie=1,ndex
nn=nv(ie)
force(2*nn)=force(2*nn)+weight*ajac*ph(1,ie)*ty
enddo
enddo
END
