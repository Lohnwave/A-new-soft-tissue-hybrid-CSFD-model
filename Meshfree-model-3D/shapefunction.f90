subroutine RPIM_ShapeFunc_2D(gpos,x,nv,phi,nx,numnode,ndex,alfc,dc,q,nRBF,mbasis)
!**************************************************************
!Compute RPIM shape function and their derivatives
!Input: gpos,x,nv,ds,alfc,dc,q,nx,numnode,ndex,mm,nRBF,nbasis
!       NRBF=1:MQ,2:EXP,3:TSP
!Output:phi
!nrbf:决定采用何形式的基函数，alfc,dc为基函数的参数；
!ndex:基函数的维数；
!mbasis:多项式函数的维数；
!numnode:节点数
!gpos:
!x:
!nv:
!ds
!alfc:
!dc:
!q:
!nx:
!From 1 to 10 of the two dimension of phi denotes
!      phi,dphidx,dphidy,dphidxx,dphidxy,dphidyy
!      dphidxxx,dphidxxy,dphidxyy,dphidyyy
!*****************************************************************
implicit real*8(a-h,o-z)
dimension gpos(nx),x(nx,numnode),nv(ndex),rk(ndex+mbasis)
dimension phi(10,ndex),xv(nx,ndex),rr(10,ndex+mbasis)
dimension a(ndex+mbasis,ndex+mbasis),g0(ndex+mbasis,ndex+mbasis),rk1(ndex+mbasis)
if(nrbf.eq.1)then
   rc=alfc*dc        !For MQ
endif
if(nrbf.eq.2)then    !For EXP
   q=alfc/dc/dc
endif

ep=1.d-20

mg=ndex+mbasis
do i=1,mg
   do j=1,mg
     g0(i,j)=0
   enddo
enddo


do i=1,ndex
   nn=nv(i)
   xv(1,i)=x(1,nn)   !xv:influence domain 中的节点数
   xv(2,i)=x(2,nn)
   xv(3,i)=x(3,nn)
enddo
!*************求整个影响域中的G0********************************
do i=1,ndex
   nn=nv(i)
   call compute_RadialBasis(x(1,nn),x(2,nn),x(3,nn),xv,rr,ndex,rc,q,nRBF,mbasis)
   do j=1,ndex
      g0(i,j)=rr(1,j)
   enddo
   if(mbasis>0)then
      g0(i,ndex+1)=1.
	  g0(i,ndex+2)=x(1,nn)
	  g0(i,ndex+3)=x(2,nn)
	  g0(ndex+1,i)=1.
	  g0(ndex+2,i)=x(1,nn)
	  g0(ndex+3,i)=x(2,nn)
    endif
enddo
!*************Solve linear equation to obtain shape function
do i=1,mg
   do j=1,mg
      a(i,j)=g0(i,j)
   enddo
enddo
!rr为输出
call compute_RadialBasis(gpos(1),gpos(2),gpos(3),xv,rr,ndex,rc,q,nRBF,mbasis)
do i=1,mg
   rk(i)=rr(1,i)
enddo
!rk既为输入，又为输出
call solverband(a,rk,mg,mg)
!  if(kwji==1)then
!   write(*,*)'Fail...'
!   pause
!     endif
    do i=1,ndex
   phi(1,i)=rk(i)
enddo

!***********************Solve linear equation to obtain dphidx
do i=1,mg
   do j=1,mg
      a(i,j)=g0(i,j)
   enddo
enddo

do i=1,mg
   rk(i)=rr(2,i)
enddo

call solverband(a,rk,mg,mg)

do i=1,ndex
   phi(2,i)=rk(i)
enddo
!***********************Solve linear equation to obtain dphidy
do i=1,mg
   do j=1,mg
      a(i,j)=g0(i,j)
   enddo
enddo

do i=1,mg
   rk(i)=rr(3,i)
enddo

call solverband(a,rk,mg,mg)

do i=1,ndex
   phi(3,i)=rk(i)
enddo
!***********************Solve linear equation to obtain dphidz
do i=1,mg
   do j=1,mg
      a(i,j)=g0(i,j)
   enddo
enddo

do i=1,mg
   rk(i)=rr(4,i)
enddo

call solverband(a,rk,mg,mg)

do i=1,ndex
   phi(4,i)=rk(i)
enddo


return
end

