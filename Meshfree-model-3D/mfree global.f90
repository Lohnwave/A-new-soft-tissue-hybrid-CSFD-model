!----------------------------------------------------------------------------
! main program--2D FORTRAN 90 CODE-MFree global weak-form methods
! Using square support domain and square background cells
! input file -- input.dat
! output file -- result.dat
! include file -- parameter.h, variable.h

implicit real*8 (a-h,o-z)
include 'parameter.h'
include 'variables.h'
open(4,file='Input175.dat')
open(2,file='result.dat',status='unknown')
!-----------测试程序运行时间-------------------------------
!Real time_begin , time_end
call CPU_TIME(time_begin)
!----------------------------------------------------------
! ************* Input data
call input(x,numd,nx,numnode,ndivx,ndivy,ndivz,ndivxq,ndivyq,ndivzq,&
nconn2,nquado,pAlf,Dmat,ALFs,numcell,numq,noCell,ncn,xc,&
npEBCnum,npEBC,pEBC,npNBCnum,npNBC,pNBC)
numgauss=nquado**3 !total number of Gauss points in a cell
! ************* Determine sizes of influence domains -- uniform nodal spacing
xspace=xlength/ndivx
yspace=ylength/ndivy
zspace=zlength/ndivz
do i=1,numnode
ds(1,i)=alfs*xspace
ds(2,i)=alfs*yspace
ds(3,i)=alfs*zspace
enddo
! ************* Coefficients of Gauss points,Weights and Jacobian for a cell
call GaussCoefficient(nquado,gauss)
do ik=1,ng
do jk=1,numgauss
gs(ik,jk)=0.
enddo
enddo
do ik=1,3*numd
force(ik)=0.
do jk=1,3*numd
ak(ik,jk)=0.
enddo
enddo
!write(2,*)'Displacements of field nodes in direction X'
!write(2,*)'Displacements of field nodes in direction Y'
!do 30 Fi=1,10,1 [5.25 -5.25 3.125]
Fx=0
Fy=0.4
Fz=0
write(*,*)'Fx=',Fx,'Fy=',Fy,'Fz=',Fz
!write(*,*)'Count No.=',Fi
! ************* Loop for background cells
do 10 ibk=1,numcell
write(*,*)'Cell No.=',ibk
! ************* Set Gauss points for this cell
call CellGaussPoints(ibk,numcell,nquado,numq,numgauss, &
xc,noCell,gauss,gs)
! ************* Loop over Gauss points to assemble discrete equations
do 20 ie=1,numgauss
gpos(1)=gs(1,ie) ! Gauss point x
gpos(2)=gs(2,ie) ! Gauss point y
gpos(3)=gs(3,ie) ! Gauss point z
weight=gs(4,ie) ! weight coefficent
ajac=gs(5,ie) ! Jacobian
! ************* Determine the support domain of Gauss point
ndex=0
call SupportDomain(numnode,nx,gpos,x,ds,ndex,nv)
do ik=1,4*ndex
do jk=1,10
ph(jk,ik)=0.
enddo
enddo
! ************* Construct RPIM shape functions for a Gauss point
call RPIM_ShapeFunc_2D(gpos,x,nv,ph,nx,numnode,ndex,&
alfc,dc,q,nRBF, mbasis)
do ik=1,3*ndex
ne(ik)=0
enddo
do ine=1,ndex
n1=3*ine-2
n2=3*ine-1
n3=3*ine
ne(n1)=3*nv(ine)-2
ne(n2)=3*nv(ine)-1
ne(n3)=3*nv(ine)
enddo
mbdb=9*ndex*ndex
do kbdb=1,mbdb
GSPk(kbdb)=0.
enddo
! ************* Compute the stiffness matrix for a Gauss point
call PointStiffnessMatrix(ndex,weight,ajac,ph,Dmat,GSPk)
nb=3*ndex
do ikk=1,nb
do jkk=1,nb
m1=ne(ikk)
m2=ne(jkk)
nbdb=(jkk-1)*nb+ikk
ak(m1,m2)=ak(m1,m2)+GSPk(nbdb)
enddo
enddo
20 continue ! end of loop for Gauss points
! ************* Implement natural BC
!in=0
!jn=0
!nn=noCell(4,ibk)
!if(xc(2,nn).eq.6) in=nn
!nn= noCell(1,ibk)
!if(xc(2,nn).eq.6) jn=nn
!if((in.ne.0).and.(jn.ne.0)) then
!call naturalBC_distributed(numnode,numq,in,jn, &
!alfs,x,xc,ds,gauss,nquado,force,Fi)
!endif
10 continue ! end of loop for cells
! ************* Boundary conditions: essential BC
write(*,*)' Boundary conditions....'
nak=2*numd
call EssentialBC(numnode,pAlf,alfs,x,ds,ak,force,npEBCnum,npEBC,pEBC)
! ************* Boundary conditions: concentrated natural BC
call NaturalBC_concentrated(x,nx,numnode,force,ds,alfs,npNBCnum,npNBC,pNBC,Fx,Fy,Fz)
nak=3*numd
b=1.d-10
! ************* Solve equation to get the solutions
!write(*,*)' Solving....',Fi/0.05
call SolverBand(ak,force,3*numnode,3*numd)
nnn=3*numd
do ik=1,nx
do jk=1, numnode
u2(ik,jk)=0.
enddo
enddo
do ik=1,numnode
jk=3*ik-2
u2(1,ik)=force(jk)
u2(2,ik)=force(jk+1)
u2(3,ik)=force(jk+2)
enddo
write(2,*)'Fx=',Fx,'Fy=',Fy,'Fz=',Fz,'top layer nodes'
do ii=1,numnode,3
write(2,52) ii,force(3*ii-2),force(3*ii-1),force(3*ii)
    enddo
write(2,*)'Fx=',Fx,'Fy=',Fy,'Fz=',Fz,'inter-layer nodes'
do ii=2,numnode,3
write(2,52) ii,force(3*ii-2),force(3*ii-1),force(3*ii)
    enddo  
write(2,*)'Fx=',Fx,'Fy=',Fy,'Fz=',Fz,'bottom layer nodes'
do ii=3,numnode,3
write(2,52) ii,force(3*ii-2),force(3*ii-1),force(3*ii)
    enddo
write(2,*)'Fx=',Fx,'Fy=',Fy,'Fz=',Fz,'all nodes'
do ii=1,numnode
write(2,52) ii,force(3*ii-2),force(3*ii-1),force(3*ii)
!-------------程序计时---------------------
call CPU_TIME(time_end)
write(2,*) 'time_begin=',time_begin
write(2,*)'time_end=',time_end
!-----------------------------------------
enddo 

52  format(1x,i5,3e20.5)

! ************* Get the final displacement
!call GetDisplacement(x,ds,u2,disp,alfs,nx,numnode,Fi)
!30 continue !end Fi
! ************* Get stress
!call GetStress(x,noCell,ds,Dmat,u2,alfs,nx,numnode,numgauss,&
!xc,gauss,nquado,ng,numq,numcell, ENORM,Stressnode)
STOP
END
