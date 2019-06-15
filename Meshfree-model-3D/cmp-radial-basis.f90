subroutine compute_RadialBasis(x,y,z,xv,rk,ndex,R,q,nRBF,mbasis)
!*********************************************************************
!Compute raidal basis after linear polynomial
!Input:x,y,xv,ndex,r,q,nRBF,mbasis
!Output:rk(10,ndex+mbasis)
!       From 1 to 10 denotes
!       r,drx,dry,drxx,dryy,drxy
!       drdxxx,drdxxy,drdxyy,drdyyy
!***********************************************************************
implicit real*8(a-h,o-z)
dimension xv(3,ndex),rk(10,ndex+mbasis)
do i=1,ndex+mbasis
   do j=1,10
      rk(j,i)=0
   enddo
enddo

do i=1,ndex
   rr2=(x-xv(1,i))**2+(y-xv(2,i))**2+(z-xv(3,i))**2
   if(nRBF.eq.1)then  !MQ
      rk(1,i)=(rr2+R**2)**q
	  rk(2,i)=2.*q*(rr2+R**2)**(q-1)*(x-xv(1,i))
	  rk(3,i)=2.*q*(rr2+R**2)**(q-1)*(y-xv(2,i))
      rk(4,i)=2.*q*(rr2+R**2)**(q-1)*(z-xv(3,i))
	  rk(5,i)=2.*q*(rr2+R**2)**(q-1)+4.*(q-1)*q*(x-xv(1,i))**2*(rr2+R**2)**(q-2)
	  rk(6,i)=4.*(q-1)*q*(x-xv(1,i))*(y-xv(2,i))*(rr2+R**2)**(q-2)
	  rk(7,i)=2.*q*(rr2+R**2)**(q-1)+4.*(q-1)*q*(y-xv(2,i))**2*(rr2+R**2)**(q-2)
    endif

	if(nRBF==2)then  !EXP
	   rk(1,i)=exp(-q*rr2)
	   rk(2,i)=-2*q*exp(-q*rr2)*(x-xv(1,i))
	   rk(3,i)=-2*q*exp(-q*rr2)*(y-xv(2,i))
	   rk(4,i)=-2*q*exp(-q*rr2)+4*q*q*(x-xv(1,i))**2*exp(-q*rr2)
	   rk(5,i)=4.*q*q*exp(-q*rr2)*(x-xv(1,i))*(y-xv(2,i))
	   rk(6,i)=-2*q*exp(-q*rr2)+4*q*q*(y-xv(2,i))**2*exp(-q*rr2)
     endif

	 if(nRBF==3)then   !TSP
	    rk(1,i)=rr2**(0.5*q)
		rk(2,i)=q*(x-xv(1,i))*rr2**(0.5*q-1)
		rk(3,i)=q*(y-xv(2,i))*rr2**(0.5*q-1)
		rk(4,i)=q*rr2**(0.5*q-1)+2.*q*(0.5*q-1)*(x-xv(1,i))**2*rr2**(0.5*q-2)
		rk(5,i)=q*(0.5*q-1)*(x-xv(1,i))*(y-xv(2,i))*rr2**(0.5*q-2)
		rk(6,i)=q*rr2**(0.5*q-1)+2.*q*(0.5*q-1)*(y-xv(2,i))**2*rr2**(0.5*q-2)
	 endif
enddo

if(mbasis>0)then
   rk(1,ndex+1)=1.
   rk(1,ndex+2)=x
   rk(1,ndex+3)=y
   rk(2,ndex+2)=1.
   rk(3,ndex+3)=1.
endif

return

end
