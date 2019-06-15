SUBROUTINE GaussSolver(n,mk,a,b,ep,kwji)
!------------------------------------------------------------------
! Stnadard Gauss elimination slover for linear equations that are
! not suitably solved by BandSolver.
!------------------------------------------------------------------
implicit real*8 (a-h,o-z)
dimension a(mk,mk),b(mk)
integer, allocatable :: m(:)
allocate (m(3*mk))
ep=1.0e-10
kwji=0
do i=1,n
m(i)=i
enddo
do 20 k=1,n
p=0.0
do 30 i=k,n
do 30 j=k,n
if(abs(a(i,j)).le.abs(p)) goto 30
p=a(i,j)
io=i
jo=j
30 continue
if(abs(p)-ep) 200,200,300
200 kwji=1
return
300 if(jo.eq.k) goto 400
do i=1,n
t=a(i,jo)
a(i,jo)=a(i,k)
a(i,k)=t
enddo
j=m(k)
m(k)=m(jo)
m(jo)=j
400 if(io.eq.k) goto 500
do j=k,n
t=a(io,j)
a(io,j)=a(k,j)
a(k,j)=t
enddo
t=b(io)
b(io)=b(k)
b(k)=t
500 p=1.0/p
in=n-1
if(k.eq.n) goto 600
do j=k,in
a(k,j+1)=a(k,j+1)*p
enddo
600 b(k)=b(k)*p
if(k.eq.n) goto 20
do i=k,in
do j=k,in
a(i+1,j+1)=a(i+1,j+1)-a(i+1,k)*a(k,j+1)
enddo
b(i+1)=b(I+1)-a(i+1,k)*b(k)
enddo
20 continue
do i1=2,n
i=n+1-i1
do j=i,in
b(i)=b(i)-a(i,j+1)*b(j+1)
enddo
enddo
do k=1,n
i=m(k)
a(1,i)=b(k)
enddo
do k=1,n
b(k)=a(1,k)
enddo
kwji=0
deallocate (m)
return
END
