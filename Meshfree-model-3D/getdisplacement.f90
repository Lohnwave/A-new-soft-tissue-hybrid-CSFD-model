SUBROUTINE GetDisplacement(x,ds,u2,disp,alfs,nx,numnode,Fi)
!------------------------------------------------------------------
! This subroutine to get the final displacements from
! displacement parameters using the MFree interpolation;
! input--numnode: total number of field nodes;
! alfs: coefficent of support support
! x(nx,numnode): coordinates of all field nodes;
! u2(2*numnode): displacement parameters;
! input and output-- disp: final displacements.
!------------------------------------------------------------------
implicit real*8 (a-h,o-z)
COMMON/rpim/ALFC,DC,Q,nRBF
common/basis/mbasis
dimension x(nx,numnode),ds(nx,numnode),gpos(nx),u2(nx,numnode)
dimension disp(3*numnode)
dimension ph(10,numnode), nv(numnode)
!write(2,*)'Displacements of field nodes'
nn=3*numnode
do i=1,nn
disp(i)=0.
enddo
ind=0
do 50 id=1,numnode
ind=ind+1
gpos(1)= x(1,id)
gpos(2)=x(2,id)
gpos(3)=x(3,id)
ndex=0
call SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)
do kph=1,ndex
do ik=1,10
ph(ik,kph)=0.
enddo
enddo
call RPIM_ShapeFunc_2D(gpos,x,nv,ph,nx,numnode,ndex,&
alfc,dc,q,nRBF, mbasis)
nc1=3*ind-2
nc2=3*ind-1
nc3=3*ind
do kk=1,ndex
m=nv(kk)
disp(nc1)=disp(nc1)+ph(1,kk)*u2(1,m)
disp(nc2)=disp(nc2)+ph(1,kk)*u2(2,m)
disp(nc3)=disp(nc3)+ph(1,kk)*u2(3,m)
enddo
50 continue
write(2,*)'Fi=',Fi
do ii=1,numnode,3
write(2,52) ii,disp(3*ii-2),disp(3*ii-1),disp(3*ii)
enddo
52 format(1x,i5,3e20.5)
!ii=1
!write(2,53)(disp(2*ii-1),disp(3*ii-1),disp(3*ii),ii=1,numnode,3)
!write(2,53)(disp(2*ii),ii=1,numnode,7)
53 format(1x,i5,363e15.5)
RETURN
END
