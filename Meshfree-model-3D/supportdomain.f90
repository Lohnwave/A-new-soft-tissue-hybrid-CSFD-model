SUBROUTINE SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)
!----------------------------------------------------------------------------
! This subroutine to determines nodes in the support domain of a Gauss point
! input--numnode: total number of field nodes;
! nx=2: for 2D problem;
! x(nx,numnode): coordinates of all field nodes;
! numgauss: number of Gauss points in a cell;
! gpos(2): x and y coordinate of a Gauss point;
! ds(nx,numnode): sizes of support domain;
! input and output-- ndex: when input ndex=0;
! when return ndex is the number of nodes in the support domain
! output--nv(ndex): No. of field nodes in the support domain
!---------------------------------------------------------------------------
implicit real*8 (a-h,o-z)
dimension gpos(nx),x(nx,numnode),ds(nx,numnode),nv(numnode)
eps=1.e-16
ndex=0
do ik=1,numnode
nv(ik)=0
enddo
do ik=1,numnode
dx=ds(1,ik)-dabs(gpos(1)-x(1,ik))
dy=ds(2,ik)-dabs(gpos(2)-x(2,ik))
dz=ds(3,ik)-dabs(gpos(3)-x(3,ik))
if((dx.ge.eps).and.(dy.ge.eps).and.(dz.ge.eps)) then
ndex=ndex+1
nv(ndex)=ik
end if
enddo
RETURN
END
