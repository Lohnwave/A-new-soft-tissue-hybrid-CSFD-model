SUBROUTINE SolverBand(ak,fp,neq,nmat)
!------------------------------------------------------------------
! Sloving linear equations; it calls BandSolver & GaussSolver
! input-ak,fp,neq,nmat
! output--fp
!------------------------------------------------------------------
implicit real*8 (a-h,o-z)
dimension ak(nmat,nEq),fp(nmat)
real(8), allocatable :: tp(:,:)
real(8), allocatable :: stfp(:,:)
allocate (tp(1:neq,1:nmat))
allocate (stfp(1:neq,1:neq))
ep=1.d-10
do i=1,nEq
do j=1,nEq
stfp(i,j)=0.
tp(i,j)=0.
enddo
enddo
do i=1,nEq
do j=1,nEq
stfp(i,j)=ak(i,j)
enddo
enddo
ni=nEq
Lp=0 ! half band width
do 20 i=1,ni
do j=ni,i,-1
if(stfp(i,j).ne.0.) then ! stfp[,] stifness matrix
if(abs(j-i).gt.Lp) Lp=abs(j-i)
go to 21
endif
enddo
21 continue
do j=1,i
if(stfp(i,j).ne.0.) then
if(abs(j-i).gt.Lp) Lp=abs(j-i)
go to 20
endif
enddo
20 continue
ilp=2*lp+1 ! band width
nm=nEq
if(ilp.lt.nEq) then
call BandSolver(stfp,fp,tp,nm,lp,ilp,nmat) ! solver for band matrix
else
call GaussSolver(nEq,nmat,ak,fp,ep,kkkk) ! standard solver
endif
deallocate (tp)
deallocate (stfp)
END

