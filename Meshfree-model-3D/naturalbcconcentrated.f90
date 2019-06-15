SUBROUTINE NaturalBC_concentrated(x,nx,numnode,af,ds,alfs,npNBCnum,npNBC,pNBC,Fx,Fy,Fz)
implicit real*8 (a-h,o-z)
dimension npNBC(4,100),pNBC(3,100)
COMMON/rpim/ALFC,DC,Q,nRBF
common/basis/mbasis
dimension af(3*numnode),x(nx,numnode),ds(nx,numnode)
dimension ph(10,numnode),gpos(3),nv(numnode)
do 10 iebc=1,npNBCnum
ie=npNBC(1,iebc)
gpos(1)=x(1,ie)
gpos(2)=x(2,ie)
gpos(3)=x(3,ie)
ndex=0
call SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)
do kph=1,3*ndex
do ik=1,10
ph(ik,kph)=0.
enddo
enddo
call RPIM_ShapeFunc_2D(gpos,x,nv,ph,nx,numnode,ndex,&
alfc,dc,q,nRBF, mbasis)
do iee=1,ndex
ie=nv(iee)
uu=pNBC(1,iebc)
af(ie*3-2)=af(ie*3-2)+ph(1,iee)*(uu-Fx)
uu=pNBC(2,iebc)
af(ie*3-1)=af(ie*3-1)+ph(1,iee)*(uu-Fy)
uu=pNBC(3,iebc)
af(ie*3)=af(ie*3)+ph(1,iee)*(uu-Fz)
enddo
10 continue
RETURN
END
